"""Check archive/extracted identity without extracting or rewriting the source."""
from pathlib import Path
from collections import Counter
import hashlib
import sys
import zipfile

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fin_r7o3_review.review import HERE,ROOT,SOURCE,load,lines,sha,save


def main():
    registry=load(HERE/'registry.json')
    archive=ROOT/(SOURCE.name+'.zip')
    if not archive.exists():
        historical=load(HERE/'provenance.json')
        assert historical['archive_manifest_match']
        assert historical['files_compared']==len(registry['source_sha256'])
        for rel,want in registry['source_sha256'].items():
            assert sha(SOURCE/rel)==want,rel
        print('ZIP absent: historical archive report preserved; all extracted inputs checked. No archive reconstructed.')
        return
    with zipfile.ZipFile(archive) as z:
        names=[n for n in z.namelist() if not n.endswith('/')]
        assert len(set(names))==len(names)
        prefixes=['',SOURCE.name+'/']
        prefix=next(p for p in prefixes if p+'MANIFEST.sha256' in names)
        assert hashlib.sha256(z.read(prefix+'MANIFEST.sha256')).hexdigest()==registry['source_manifest_sha256']
        for rel,want in registry['source_sha256'].items():
            assert hashlib.sha256(z.read(prefix+rel)).hexdigest()==want,rel
    certs=lines(SOURCE/'certificates/active_leaf_certificates.jsonl')
    declared=Counter(c['checker_source_sha256'] for c in certs)
    actual=sha(SOURCE/'src/fixed_witness_checker.py')
    save('provenance.json',dict(archive_sha256=sha(archive),archive_path=archive.name,
          archive_manifest_match=True,files_compared=len(registry['source_sha256']),
          shipped_checker_sha256=actual,declared_checker_sha256=dict(declared),
          declared_checker_matches_shipped=all(v==actual for v in declared),
          resolution='Historical checker hash not promoted as current provenance. Mathematical acceptance uses fresh replay of the inspected shipped checker with its actual hash.'))
    print('ZIP/extracted manifest and',len(registry['source_sha256']),'files match; historical checker metadata recorded')


if __name__=='__main__':
    if not __debug__:raise RuntimeError('Assertions must be enabled')
    main()
