# calling R from python
import subprocess
import sys

command = 'Rscript'

if __name__ == '__main__':
    path2script = "./simple.R"

    cmd = [command, path2script] + sys.argv[1:]
    print("cmd is %s" %cmd)

    x = subprocess.check_output(cmd, universal_newlines=True)
    print("cmd output is: ")
    print("%s" %str(x))

