import subprocess

def AeroElastic():
    #subprocess.call(['sh', './execute.sh'])
    result = subprocess.Popen([ "./execute.sh"])
    result.wait()
    #text = result.communicate()[0]
    returncode = result.returncode
    return returncode
