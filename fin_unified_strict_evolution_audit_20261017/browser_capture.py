"""Run inside browser-harness; public page provenance, not a theorem checker."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlparse

root = Path.cwd()
info = page_info()
assert urlparse(info['url']).hostname == 'arxiv.org'
assert all(c.isalnum() or c in '-_' for c in CAPTURE_LABEL)
body = js('document.body.innerText')
assert isinstance(body, str) and len(body) > 200
cache = root / '.firecrawl/fin_unified_evolution_20261017'
cache.mkdir(parents=True, exist_ok=True)
path = cache / (CAPTURE_LABEL + '.txt')
path.write_text(body)
entry = {'label': CAPTURE_LABEL, 'observed_utc': datetime.now(timezone.utc).isoformat(),
         'url': info['url'], 'title': info['title'],
         'method': 'actual existing local browser-harness, recording off',
         'text_sha256': hashlib.sha256(body.encode()).hexdigest(),
         'text_bytes': len(body.encode()), 'optional_cache': str(path.relative_to(root)),
         'scope': 'public DOM access and discovery, not full PDF proof verification'}
with (root/'fin_unified_strict_evolution_audit_20261017/BROWSER_READ_LOG.jsonl').open('a') as f:
    f.write(json.dumps(entry, ensure_ascii=False)+'\n')
print(json.dumps(entry, ensure_ascii=False))
