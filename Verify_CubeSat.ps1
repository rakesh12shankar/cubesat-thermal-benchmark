$ErrorActionPreference='Stop'
$interop='C:\Program Files\SOLIDWORKS Corp\SOLIDWORKS\SolidWorks.Interop.sldworks.dll'
Add-Type -Path $interop
Add-Type -ReferencedAssemblies $interop -TypeDefinition @'
using System;
using System.IO;
using SolidWorks.Interop.sldworks;
public class CubeSatVerify {
 public static string Run(string path){
  ISldWorks app=(ISldWorks)Activator.CreateInstance(Type.GetTypeFromProgID("SldWorks.Application"));
  int error=0,warning=0;
  IModelDoc2 doc=(IModelDoc2)app.OpenDoc6(path,1,1,"",ref error,ref warning);
  if(doc==null||error!=0)throw new Exception("Cannot reopen saved CAD: "+error);
  object[] bodies=(object[])((IPartDoc)doc).GetBodies2(0,false);
  if(bodies.Length!=27)throw new Exception("Incorrect saved body count");
  double total=0;int pairs=0;
  for(int i=0;i<bodies.Length;i++){
   IBody2 body=(IBody2)bodies[i];total+=((double[])body.GetMassProperties(1))[3];
   for(int j=i+1;j<bodies.Length;j++){
    int err;object raw=((IBody2)body.Copy()).Operations2(15901,((IBody2)bodies[j]).Copy(),out err);
    object[] intersection=raw as object[];
    if(intersection!=null)foreach(IBody2 overlap in intersection){
     double volume=((double[])overlap.GetMassProperties(1))[3];
     if(volume>1e-12)throw new Exception("Volumetric overlap between bodies "+i+" and "+j);
    }
    pairs++;
   }
  }
  if(Math.Abs(total-0.000247936)>1e-12)throw new Exception("Saved volume mismatch: "+total);
  return "Saved native CAD reopened; 27 solids; "+pairs+" intersection checks; total solid volume="+total.ToString("R",System.Globalization.CultureInfo.InvariantCulture)+" m3.";
 }
}
'@
$result=[CubeSatVerify]::Run((Join-Path $PSScriptRoot 'data/cad/CubeSat_2021_nominal.SLDPRT'))
$result | Set-Content -LiteralPath (Join-Path $PSScriptRoot 'data/cad/CAD_Verification.txt')
$result
