"""Bundle preserved records after completion; never alter the source study."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import zipfile


def sha(data):
    return hashlib.sha256(data).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    study, out = args.study.resolve(), args.out.resolve()
    if study.drive.lower() != 'd:' or out.drive.lower() != 'd:':
        raise ValueError('Keep study and archive output on D')
    if out.is_relative_to(study):
        raise ValueError('Keep the archive outside its source study')
    repository = str(Path(__file__).resolve().parents[2])
    if repository not in sys.path:
        sys.path.insert(0, repository)
    from d3d_optimization.campaign_state import process_lock
    with process_lock(study/'campaign.lock'), process_lock(study/'posttrial.lock'):
        build(study, out)


def build(study, out):
    for phase in ('reference', 'uncertainty', 'screen-rejections'):
        if not (study / phase / 'summary.json').is_file():
            raise ValueError('Wait for the physical post-trial stage to finish: ' + phase)
    frozen = json.loads((study / 'frozen.json').read_text(encoding='utf-8'))
    amendment = json.loads((study / 'amendment.json').read_text(encoding='utf-8'))
    prior = Path(amendment['prior_study']).resolve()
    blocked = {'cache', 'posttrial-runtime', '__pycache__', '.matplotlib'}
    members, exclusions = {}, []
    for root in (study, prior):
        for path in sorted(root.rglob('*')):
            if not path.is_file():
                continue
            relative = path.relative_to(root)
            if blocked.intersection(relative.parts) or path.suffix == '.lock':
                exclusions.append(str(path))
                continue
            name = 'studies/' + root.name + '/' + relative.as_posix()
            if name in members:
                raise ValueError('Duplicate archive member: ' + name)
            members[name] = path
    source = Path(frozen['source_study'])
    for name in ('selection.json', 'downloads.json', 'preflight.json', 'protocol.ko.md'):
        members['studies/' + source.name + '/' + name] = source / name
    entries = {name: dict(bytes=path.stat().st_size, sha256=sha(path.read_bytes()))
               for name, path in sorted(members.items())}
    manifest = dict(schema='nf-anytime-evidence-bundle.v1', study=str(study),
                    entries=entries, member_count=len(entries),
                    uncompressed_bytes=sum(row['bytes'] for row in entries.values()),
                    excluded_runtime_files=exclusions,
                    external_input_scope='Original public parquet files and prepared case input roots remain external. Frozen input hashes and download URLs are retained; see reproduce.ko.md.',
                    restoration_scope='Inspect/extract into a fresh directory. This archive does not install dependencies or replace an existing study.')
    payload = (json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False)+'\n').encode('utf-8')
    out.mkdir(parents=True, exist_ok=True)
    archive = out / 'evidence.zip'
    manifest_path = out / 'evidence-manifest.json'
    checksum_path = out / 'evidence.zip.sha256'
    if any(p.exists() for p in (archive, manifest_path, checksum_path)):
        raise FileExistsError('Choose a fresh output directory; existing evidence is never overwritten')
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as bundle:
        for name, path in sorted(members.items()):
            data = path.read_bytes()
            if len(data) != entries[name]['bytes'] or sha(data) != entries[name]['sha256']:
                raise ValueError('Source changed while bundling: ' + str(path))
            info = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            bundle.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
        info = zipfile.ZipInfo('evidence-manifest.json', date_time=(1980, 1, 1, 0, 0, 0))
        info.compress_type = zipfile.ZIP_DEFLATED
        info.external_attr = 0o100644 << 16
        bundle.writestr(info, payload, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(archive) as bundle:
        if set(bundle.namelist()) != set(entries) | {'evidence-manifest.json'}:
            raise ValueError('Unexpected archive members')
        for name, row in entries.items():
            data = bundle.read(name)
            if len(data) != row['bytes'] or sha(data) != row['sha256']:
                raise ValueError('Archive round-trip differs: ' + name)
        if bundle.read('evidence-manifest.json') != payload:
            raise ValueError('Archived manifest differs')
    with manifest_path.open('xb') as stream:
        stream.write(payload)
    with checksum_path.open('x', encoding='utf-8') as stream:
        stream.write(sha(archive.read_bytes())+'  evidence.zip\n')
    print(json.dumps(dict(members=len(entries), archive_bytes=archive.stat().st_size,
                          uncompressed_bytes=manifest['uncompressed_bytes'],
                          every_archived_member_sha256_verified=True)))


if __name__ == '__main__':
    main()
