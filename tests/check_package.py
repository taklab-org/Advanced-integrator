"""Source-only release audit; uses Python standard library only."""
from pathlib import Path
import hashlib,json,re
root=Path(__file__).resolve().parents[1]
files=[p for p in root.rglob('*') if p.is_file() and '.git' not in p.parts and 'outputs' not in p.relative_to(root).parts]
for p in files:
 assert p.suffix.lower() not in {'.mat','.log','.dll','.mexw64','.so','.dylib','.zip'},p
 if p.suffix=='.m':
  s=p.read_text(encoding='utf-8-sig')
  assert not re.search(r"(?:[A-Z]:\\|/Users/|/home/)\w",s),p
s=(root/'cases/num7/proof_GE_2DOK_core.m').read_text(encoding='utf-8-sig')
assert 'tau = 8;' in s and 'K_F_list = 128;' in s and 'm_a=[12,12];h_auto=0.1' in s
assert '8.625' not in s and 'K_F_list = 138' not in s
assert len(list(root.glob('cases/num*/initial_num*.m')))==2
manifest=[]
for p in sorted(files):
 if p.relative_to(root).as_posix()=='docs/release-files.json':continue
 manifest.append({'path':p.relative_to(root).as_posix(),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(root/'docs/release-files.json').write_text(json.dumps(manifest,indent=2)+'\n',encoding='utf-8',newline='\n')
print(f'PASS: {len(manifest)} source/document files; no binary computation data or host paths.')
