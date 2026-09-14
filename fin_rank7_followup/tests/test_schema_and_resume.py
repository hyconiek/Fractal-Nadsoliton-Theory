import json,sys,tempfile,unittest,subprocess,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from src.schema import validate_claim,validate_interval,validate_certificate
from src.resumable_cover import run

class TestSchemas(unittest.TestCase):
    def test_claim_missing_domain_rejected(self):
        with self.assertRaises(ValueError): validate_claim({'id':'x','statement':'x','status':'UNRESOLVED','dimension':'1','finite_or_limit':'finite','source':'x','evidence_level':'x'})
    def test_float_exact_endpoint_rejected(self):
        with self.assertRaises(ValueError): validate_interval({'lo':0.1,'hi':0.2},certified=True)
    def test_rational_interval_accepts(self):
        self.assertTrue(validate_interval({'lo':'1/10','hi':{'num':1,'den':5}},certified=True))
    def test_global_pass_with_unresolved_rejected(self):
        c={'id':'C','claim_id':'X','domain':'D','quantifiers':'all','assumptions':['a'],'proof_type':'cover','inputs':['i'],'conclusion':'x','global_pass':True,'unresolved_leaves':[1]}
        with self.assertRaises(ValueError): validate_certificate(c)

class TestResume(unittest.TestCase):
    def test_resume_matches_uninterrupted(self):
        with tempfile.TemporaryDirectory() as d:
            a=Path(d)/'a.json'; b=Path(d)/'b.json'
            full=run(a)
            part=run(b,stop_after=17)
            self.assertFalse(part['complete'])
            resumed=run(b)
            self.assertTrue(resumed['complete'])
            self.assertEqual(full['leaves'],resumed['leaves'])
            self.assertEqual(full['processed'],resumed['processed'])

class TestMutation(unittest.TestCase):
    def test_input_hash_mutation_rejected(self):
        # Work in copied package skeleton to avoid touching real inputs.
        with tempfile.TemporaryDirectory() as d:
            dst=Path(d)/'pkg'; shutil.copytree(ROOT,dst)
            cp=subprocess.run([sys.executable,str(dst/'verify.py')],cwd=dst,text=True,capture_output=True)
            self.assertEqual(cp.returncode,0,cp.stdout+cp.stderr)
            target=next((dst/'inputs').glob('FIN_Post_Handoff_Research_Master_Plan_EN.md'))
            target.chmod(0o644); target.write_text(target.read_text()+'\nMUTATION\n')
            cp=subprocess.run([sys.executable,str(dst/'verify.py')],cwd=dst,text=True,capture_output=True)
            self.assertNotEqual(cp.returncode,0)
            self.assertIn('hash mismatch',cp.stdout)

if __name__=='__main__':unittest.main()
