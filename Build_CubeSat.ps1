$ErrorActionPreference='Stop'
$interop='C:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\SolidWorks.Interop.sldworks.dll'
Add-Type -Path $interop
Add-Type -ReferencedAssemblies $interop -TypeDefinition @'
using System;
using System.IO;
using System.Collections.Generic;
using System.Globalization;
using SolidWorks.Interop.sldworks;
public class CubeSatCad {
 static IModeler modeler;
 static IBody2 Box(double x,double y,double z,double a,double b,double c){
  IBody2 body=(IBody2)modeler.CreateBodyFromBox(new double[]{(x+a/2)/1000,(y+b/2)/1000,z/1000,0,0,1,a/1000,b/1000,c/1000});
  if(body==null)throw new Exception("Box failed"); return body;
 }
 static IBody2 Sub(IBody2 a,IBody2 b){
  int error;object[] result=(object[])a.Operations2(15902,b.Copy(),out error);
  if(error!=0||result==null||result.Length!=1)throw new Exception("Subtract failed: "+error);return (IBody2)result[0];
 }
 public static string Run(string directory){
  Directory.CreateDirectory(directory);string stem=Path.Combine(directory,"CubeSat_2021_nominal");
  if(File.Exists(stem+".SLDPRT"))throw new Exception("CAD exists; refusing overwrite");
  ISldWorks app=(ISldWorks)Activator.CreateInstance(Type.GetTypeFromProgID("SldWorks.Application"));
  modeler=(IModeler)app.GetModeler();
  var bodies=new List<IBody2>();var names=new List<string>();
  bodies.Add(Box(0,0,0,100,100,2));names.Add("Panel_5_bottom");
  bodies.Add(Box(0,0,98,100,100,2));names.Add("Panel_6_top");
  bodies.Add(Box(0,0,2,2,100,96));names.Add("Panel_1_minusX");
  bodies.Add(Box(98,0,2,2,100,96));names.Add("Panel_2_plusX_nadir");
  bodies.Add(Box(2,0,2,96,2,96));names.Add("Panel_3_minusY");
  bodies.Add(Box(2,98,2,96,2,96));names.Add("Panel_4_plusY");
  var metal=new List<IBody2>();int index=0;
  foreach(double x in new double[]{2,93})foreach(double y in new double[]{2,93}){
   var v=Box(x,y,2,5,5,96);metal.Add(v);bodies.Add(v);names.Add("Frame_vertical_"+(++index));
  }
  foreach(double z in new double[]{2,93}){
   foreach(double y in new double[]{2,93}){var v=Box(7,y,z,86,5,5);metal.Add(v);bodies.Add(v);names.Add("Frame_X_"+(++index));}
   foreach(double x in new double[]{2,93}){var v=Box(x,7,z,5,86,5);metal.Add(v);bodies.Add(v);names.Add("Frame_Y_"+(++index));}
  }
  index=0;
  foreach(double x in new double[]{10,85})foreach(double y in new double[]{10,85}){
   var v=Box(x,y,2,5,5,96);metal.Add(v);bodies.Add(v);names.Add("Bolt_"+(++index));
  }
  index=0;
  foreach(double z in new double[]{20,37,61,78}){
   var v=Box(5,5,z,90,90,2);
   // Cut the four corner-frame intrusions and bolt holes; one board remains.
   foreach(double x in new double[]{3,93})foreach(double y in new double[]{3,93})v=Sub(v,Box(x,y,z-1,4,4,4));
   foreach(double x in new double[]{10,85})foreach(double y in new double[]{10,85})v=Sub(v,Box(x,y,z-1,5,5,4));
   bodies.Add(v);names.Add("PCB_"+(++index));
  }
  bodies.Add(Box(20,20,39,60,60,9));names.Add("Battery");
  IModelDoc2 doc=(IModelDoc2)app.NewDocument(@"C:\ProgramData\SOLIDWORKS\SOLIDWORKS 2019\templates\Part.prtdot",0,0,0);
  if(doc==null)throw new Exception("NewDocument failed");IPartDoc part=(IPartDoc)doc;
  var rows=new List<string>();rows.Add("body,volume_m3");
  for(int i=0;i<bodies.Count;i++){
   IFaultEntity fault=(IFaultEntity)bodies[i].Check3;if(fault!=null&&fault.Count>0)throw new Exception("Invalid body: "+names[i]);
   double volume=((double[])bodies[i].GetMassProperties(1))[3];
   IFeature feature=(IFeature)part.CreateFeatureFromBody3(bodies[i].Copy(),false,1);
   if(feature==null)throw new Exception("Feature failed");feature.Name=names[i];
   rows.Add(names[i]+","+volume.ToString("R",CultureInfo.InvariantCulture));
  }
  object[] saved=(object[])part.GetBodies2(0,false);if(saved.Length!=27)throw new Exception("Expected 27 solids");
  doc.ForceRebuild3(false);doc.ShowNamedView2("*Isometric",7);doc.ViewZoomtofit2();
  foreach(string extension in new string[]{"SLDPRT","x_t","STEP"}){
   int error=0,warning=0;app.ActivateDoc3(doc.GetTitle(),false,0,ref error);error=0;
   if(!doc.Extension.SaveAs(stem+"."+extension,0,1,null,ref error,ref warning)||error!=0)throw new Exception("Save failed: "+error);
  }
  File.WriteAllLines(Path.Combine(directory,"CAD_Volumes.csv"),rows);return "Saved 27 valid solids in SLDPRT, Parasolid, and STEP formats.";
 }
}
'@
[CubeSatCad]::Run((Join-Path $PSScriptRoot 'data/cad'))

