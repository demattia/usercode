import math

def deltaPhi(phi1, phi2):
    result = phi1 - phi2;
    if result > math.pi: result -= 2*math.pi
    if result <= -math.pi: result += 2*math.pi
    return result

def deltaR(v1, v2):
    return math.sqrt(deltaPhi(v1.phi(), v2.phi())**2 + (v1.eta()-v2.eta())**2)

def deltaR(phi1, eta1, phi2, eta2):
    return math.sqrt(deltaPhi(phi1, phi2)**2 + (eta1-eta2)**2)

# This is for intercepting the output of ROOT
# In a cell, put %%rootprint so that the output that would normally be
# sent directly to the stdout will instead be displayed in the cell.
# It must be the first element in the cell.
import tempfile
import ROOT
from IPython.core.magic import (Magics, magics_class, cell_magic)

@magics_class
class RootMagics(Magics):
    """Magics related to Root.

    %%rootprint  - Capture Root stdout output and show in result cell
    """

    def __init__(self, shell):
        super(RootMagics, self).__init__(shell)

    @cell_magic
    def rootprint(self, line, cell):
        """Capture Root stdout output and print in ipython notebook."""

        with tempfile.NamedTemporaryFile() as tmpFile:

            ROOT.gSystem.RedirectOutput(tmpFile.name, "w")
            ns = {}
            exec cell in self.shell.user_ns, ns
            ROOT.gROOT.ProcessLine("gSystem->RedirectOutput(0);")
            print tmpFile.read()

# Register
ip = get_ipython() 
ip.register_magics(RootMagics)
