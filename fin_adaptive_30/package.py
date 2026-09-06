"""Create a checked portable archive; optionally verify a clean extraction."""
import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent


def main():
    manifest=HERE/'SHA256SUMS.txt'
    names=[]
    for row in manifest.read_text().splitlines():
        expected,name=row.split('  ',1)
        if Path(name).name != name:
            raise ValueError('Only local package files may be archived')
        actual=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        if actual != expected:
            raise ValueError(f'Stale manifest for {name}; run build_artifacts.py')
        names.append(name)
    target=HERE.parent/'FIN_ST8549_ST8578_Thirty_Adaptive_Studies_Package.zip'
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as archive:
        for name in names+['SHA256SUMS.txt']:
            archive.write(HERE/name,arcname='fin_adaptive_30/'+name)
    digest=hashlib.sha256(target.read_bytes()).hexdigest()
    target.with_suffix('.sha256').write_text(digest+'  '+target.name+'\n')
    with zipfile.ZipFile(target) as archive:
        if archive.testzip() is not None:
            raise RuntimeError('ZIP integrity failure')
        if '--verify' in sys.argv:
            with tempfile.TemporaryDirectory(prefix='fin30_package_') as tmp:
                archive.extractall(tmp)  # paths were generated and checked above
                work=Path(tmp)/'fin_adaptive_30'
                subprocess.run(['sha256sum','-c','SHA256SUMS.txt'],cwd=work,check=True)
                subprocess.run([sys.executable,'verify.py'],cwd=work,check=True)
    print(target)
    print('SHA256 '+digest)


if __name__=='__main__':
    main()
