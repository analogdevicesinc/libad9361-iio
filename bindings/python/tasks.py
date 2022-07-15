from invoke import task

@task
def bumpversion_test(c):
    """Bump version to {current-version}.dev.{date}
    Used for marking development releases for test-pypi
    """
    import fileinput
    import time

    for line in fileinput.input("setup.py", inplace=True):
        if line.find("version") > -1:
            l = line[len("version="):].strip()[:-1].strip("\"=version")[:].split(".")
            major = (l[0])
            minor = (l[1])
            seconds = int(time.time())
            line = '        version="{}.{}.{}",\n'.format(
                major, minor, seconds
            )
            ver_string = "v{}.{}".format(major, minor, seconds)
        print(line, end="")

    print(f"Version bumped to {ver_string}")
