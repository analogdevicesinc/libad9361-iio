from invoke import task

@task
def bumpversion_test(c):
    """Bump version to {current-version}.{timestamp}
    Used for marking development releases for test-pypi
    """
    import fileinput
    import time
    import re
    import os

    pyproject_file = os.path.join("pyproject.toml")

    if not os.path.exists(pyproject_file):
        print(f"Error: {pyproject_file} not found. Run cmake first.")
        return

    ver_string = None

    for line in fileinput.input(pyproject_file, inplace=True):
        if line.startswith("version"):
            match = re.search(r'version\s*=\s*"([^"]+)"', line)
            if match:
                version = match.group(1)
                parts = version.split(".")
                major = parts[0]
                minor = parts[1] if len(parts) > 1 else "0"
                seconds = int(time.time())
                new_version = f"{major}.{minor}.{seconds}"
                line = f'version = "{new_version}"\n'
                ver_string = f"v{new_version}"
        print(line, end="")

    if ver_string:
        print(f"Version bumped to {ver_string}")
