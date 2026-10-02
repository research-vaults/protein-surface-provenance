#!/usr/bin/env python3
"""Offline cache-path tests; no cache creation or downloads."""
import json
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPRODUCTION = os.path.join(ROOT, "reproduction")
CAPACITY = os.path.join(REPRODUCTION, "reproduce_channel_capacity.py")
FROZEN = os.path.join(REPRODUCTION, "reproduce_frozen_contrast.py")
DEFAULT_CACHE = os.path.join(REPRODUCTION, ".cache")
checks = 0


def resolved(script, arguments=(), cache_env=None):
    env = os.environ.copy()
    if cache_env is None:
        env.pop("SURFACE_PROVENANCE_CACHE", None)
    else:
        env["SURFACE_PROVENANCE_CACHE"] = cache_env
    result = subprocess.run(
        [sys.executable, script, "--show-paths", *arguments],
        check=True,
        capture_output=True,
        text=True,
        env=env,
    )
    return json.loads(result.stdout)


def check(label, condition):
    global checks
    if not condition:
        raise AssertionError(label)
    checks += 1
    print(f"PASS  {label}")


with tempfile.TemporaryDirectory() as directory:
    env_cache = os.path.join(directory, "environment cache")
    explicit_data = os.path.join(directory, "explicit data")
    explicit_cache = os.path.join(directory, "explicit cache")

    got = resolved(CAPACITY)
    check("capacity default", got["data_dir"] == os.path.join(DEFAULT_CACHE, "data"))

    got = resolved(CAPACITY, cache_env=env_cache)
    check(
        "capacity environment",
        got["data_dir"] == os.path.join(env_cache, "data")
        and not os.path.exists(env_cache),
    )

    got = resolved(CAPACITY, ("--data-dir", explicit_data), env_cache)
    check(
        "capacity explicit override",
        got["data_dir"] == explicit_data and not os.path.exists(explicit_data),
    )

    got = resolved(FROZEN)
    check("frozen default root", got["cache_dir"] == DEFAULT_CACHE)
    check("frozen default data", got["data_dir"] == os.path.join(DEFAULT_CACHE, "data"))

    got = resolved(FROZEN, cache_env=env_cache)
    check(
        "frozen environment root",
        got["cache_dir"] == env_cache and not os.path.exists(env_cache),
    )
    check("frozen environment data", got["data_dir"] == os.path.join(env_cache, "data"))

    got = resolved(FROZEN, ("--cache-dir", explicit_cache), env_cache)
    check(
        "frozen explicit root",
        got["cache_dir"] == explicit_cache and not os.path.exists(explicit_cache),
    )
    check("frozen explicit data", got["data_dir"] == os.path.join(explicit_cache, "data"))

print(f"ALL {checks} CACHE-PATH TESTS PASSED")
