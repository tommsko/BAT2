import os
import hashlib
import requests
from pathlib import Path

API_BASE = "https://bioid.dyn.cloud.e-infra.cz:443"
LOCAL_DIR = "."


def local_hash(path, block_size=65536):
    hasher = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(block_size), b""):
            hasher.update(block)
    return hasher.hexdigest()


def ensure_dir(path):
    path = os.path.dirname(path)
    if not path or path == '.':
        return
    os.makedirs(path, exist_ok=True)


def sync():
    print("[SYNC] Fetching database list from API...")
    resp = requests.get(f"{API_BASE}/api/LTS", timeout=30)
    resp.raise_for_status()
    remote_files = resp.json()

    for rf in remote_files:
        rel_path = rf["path"]
        remote_hash = rf["hash"]

        local_path = Path(LOCAL_DIR) / rel_path
        if not local_path.exists():
            print(f"... Downloading new file: {rel_path}")
        else:
            local_hash_value = local_hash(local_path)
            if local_hash_value == remote_hash:
                print(f"... File up to date: {rel_path}")
                continue
            print(f"... Updating changed file: {rel_path}")

        ensure_dir(local_path)
        r = requests.get(f"{API_BASE}/api/LTS/download/{rel_path.replace('/', '_').replace('.', '|')}", stream=True)
        try:
            r.raise_for_status()
            with open(local_path, "wb") as f:
                for chunk in r.iter_content(chunk_size=8192):
                    f.write(chunk)
        except Exception as exc:
            print(f"... FAILED: {exc}")

    print("[SYNC] Sync complete.")
