#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
GEN = Path(__file__).resolve().parent / "generated"
ROOT = Path(__file__).resolve().parent.parent

def _read(p: Path) -> dict:
    return json.loads(p.read_text(encoding='utf-8'))

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument('--p1310', type=Path, default=GEN / 'p1310_qw2191_r26_nb1_postclosure_independent_replay_audit_summary.json')
    ap.add_argument('--release', type=Path, default=ROOT / 'RELEASE_8_NB_CLOSURE_AUDIT_EN_PL.md')
    ap.add_argument('--out', type=Path, default=GEN / 'p1311_qw2191_r27_nb1_archival_and_external_disclosure_packet_summary.json')
    a = ap.parse_args()
    p1310 = _read(a.p1310)
    if p1310.get('next_priority') != 'R27_NB1_ARCHIVAL_AND_EXTERNAL_DISCLOSURE_PACKET':
        raise SystemExit('P1311 requires next_priority=R27_NB1_ARCHIVAL_AND_EXTERNAL_DISCLOSURE_PACKET from P1310.')
    if p1310.get('r26_independent_replay_audit', {}).get('status') != 'AUDIT_COMPLETE':
        raise SystemExit('P1311 requires AUDIT_COMPLETE from P1310.')
    if not a.release.exists():
        raise SystemExit(f'Missing successor release file: {a.release}')
    out = {
        'packet': 'P1311', 'as_of': '2026-05-11', 'lane': 'NB1_NONBRIDGE_TRACK',
        'input': {'p1310': str(a.p1310), 'release': str(a.release)},
        'disclosure': {'status': 'READY', 'audited': True, 'scope': 'NB_ONLY'},
        'next_priority': 'R28_NB1_LONG_HORIZON_MONITORING_PACKET'
    }
    a.out.write_text(json.dumps(out, indent=2, sort_keys=True)+'\n', encoding='utf-8')
    print(f"[P1311] wrote {a.out}; disclosure={out['disclosure']['status']}")

if __name__ == '__main__':
    main()
