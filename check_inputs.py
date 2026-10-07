"""Fast checks of geometry, units, load geometry and input safety; no solver."""
import importlib.util
from pathlib import Path
import tempfile
import unittest
import numpy as np
spec=importlib.util.spec_from_file_location('model',Path(__file__).with_name('prepare_model.py'))
model=importlib.util.module_from_spec(spec);spec.loader.exec_module(model)

class InputChecks(unittest.TestCase):
    def test_geometry_integrals(self):
        xyz,elems,parts,faces=model.make_mesh(15)
        self.assertAlmostEqual(sum(e[4] for e in elems),.000247936,places=12)
        self.assertAlmostEqual(sum(e[4]*model.PROPS[e[2]][0] for e in elems),.57848424,places=10)
        for side in range(1,7):self.assertAlmostEqual(sum(f[3] for f in faces if f[2]==side),.01,places=12)
        self.assertTrue(np.all((xyz>=0)&(xyz<=.1)))
        self.assertEqual(len(parts),13)

    def test_earth_view_geometry(self):
        flux,earth=model.flux_history()
        # Nadir disk factor = sin(theta)^2; opposite face cannot see Earth.
        self.assertAlmostEqual(earth[1],(6371/(6371+431))**2,places=10)
        self.assertEqual(earth[0],0)
        np.testing.assert_allclose(earth[2:],earth[2],atol=1e-12)
        self.assertTrue(np.all(flux[:,6:8]==0)) # ±Z sees no direct Sun in beta=0 orbit.
        np.testing.assert_allclose(flux[0,2:],flux[-1,2:],atol=1e-10)

    def test_eclipse_blocks_direct_sun(self):
        flux,_=model.flux_history()
        eclipse=flux[:,1].astype(bool)
        self.assertTrue(eclipse.any())
        self.assertTrue(np.all(flux[eclipse,2:8]==0))
        self.assertTrue(np.all(flux[:,14:20]>=0))
        np.testing.assert_allclose(flux[:,14:20],np.broadcast_to(flux[0,14:20],flux[:,14:20].shape))

    def test_invalid_inputs_and_overwrite_protection(self):
        for h,eps in [(0,.5),(-1,.5),(15,1.1),(15,-.1)]:
            with self.assertRaises(ValueError):model.prepare(h,eps,100,'cooldown')
        old=model.ROOT
        with tempfile.TemporaryDirectory() as folder:
            model.ROOT=Path(folder);(model.ROOT/'data/ansys').mkdir(parents=True)
            run=model.ROOT/'data/ansys/cooldown_h15_e0.5';run.mkdir();(run/'keep.txt').write_text('retained')
            try:
                with self.assertRaises(FileExistsError):model.prepare(15,.5,100,'cooldown')
                self.assertEqual((run/'keep.txt').read_text(),'retained')
            finally:model.ROOT=old

if __name__=='__main__':unittest.main()
