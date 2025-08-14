#!/usr/bin/env python3
import os
import subprocess
import stat

def is_executable_file(path):
    """Return True if `path` is a regular file and is marked executable."""
    try:
        st = os.stat(path)
    except OSError:
        return False
    return stat.S_ISREG(st.st_mode) and bool(st.st_mode & (stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH))

def run_all_in(dir_name):
    base = os.path.join(os.getcwd(), dir_name)
    if not os.path.isdir(base):
        print(f"Directory not found: {base}")
        return
    print(f"\n=== Running in {dir_name}/ ===")
    for entry in sorted(os.listdir(base)):
        full = os.path.join(base, entry)
        # skip .cpp files, directories, etc.
        if entry.endswith('.cpp') or entry.startswith('.') or not is_executable_file(full):
            continue
        print(f"---> {entry}")
        # run with cwd so that any relative paths in the program still work
        subprocess.run([f"./{entry}"], cwd=base, check=True)

def main():
    for sub in ("INCR", "NAIVE"):
        run_all_in(sub)

if __name__ == "__main__":
    main()
