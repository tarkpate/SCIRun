import unittest
import subprocess
import os

def getSourceRoot():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
    
def getNetworkFolder():
    return os.path.join(getSourceRoot(), "ExampleNets", "regression")

def getSCIRunExecutable():
    return os.path.join(os.path.dirname(getSourceRoot()), "bin", "SCIRun", "Release", "SCIRun.exe")

    
class TestSCIRunAddModule(unittest.TestCase):

    def test_scirun_module_ids(self):
        self.assertTrue(True)
        self.assertEqual(getSCIRunExecutable(), r"C:\Dev\SCIRun\bin\SCIRun\Release\SCIRun.exe")
        self.assertEqual(getSourceRoot(), r"C:\Dev\SCIRun\src")
        self.assertEqual(getNetworkFolder(), r"C:\Dev\SCIRun\src\ExampleNets\regression")
        subprocess.check_output([r"C:\_\SCIRun\Release\SCIRun.exe", os.path.join(getNetworkFolder(), "basicViz.srn5"), "-E", "--no_splash"])
        #subprocess.check_output([r"C:\_\SCIRun\Release\SCIRun.exe", r"C:\Dev\SCIRun\src\ExampleNets\regression\setuptdcs_colin27_patchelc.srn5", "-E", "--no_splash"])
        

if __name__ == '__main__':
    unittest.main()
