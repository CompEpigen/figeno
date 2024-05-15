import os
from figeno import figeno_make
import matplotlib.testing.compare

os.chdir("test_data")
def test_nanopore():
    figeno_make(config_file="GDM1_config.json",config={"output":{"file":"GDM1_figure2.svg"}})
    assert matplotlib.testing.compare.compare_images("GDM1_figure.svg","GDM1_figure2.svg",tol=0.1) is None

def test_nanopore_splitreads():
    figeno_make(config_file="GDM1_splitread_config.json",config={"output":{"file":"GDM1_splitread_figure2.svg"}})
    assert matplotlib.testing.compare.compare_images("GDM1_splitread_figure.svg","GDM1_splitread_figure2.svg",tol=0.1) is None

def test_HiC():
    figeno_make(config_file="LNCaP_config.json",config={"output":{"file":"LNCaP_figure2.svg"}})
    assert matplotlib.testing.compare.compare_images("LNCaP_figure.svg","LNCaP_figure2.svg",tol=0.1) is None

def test_wgs_circos():
    figeno_make(config_file="THP1_circos_config.json",config={"output":{"file":"THP1_circos_figure2.svg"}})
    assert matplotlib.testing.compare.compare_images("THP1_circos_figure.svg","THP1_circos_figure2.svg",tol=0.1) is None

def test_wgs_symmetrical():
    figeno_make(config_file="THP1_symmetrical_config.json",config={"output":{"file":"THP1_symmetrical_figure2.png"}})
    assert matplotlib.testing.compare.compare_images("THP1_symmetrical_figure.png","THP1_symmetrical_figure2.png",tol=0.1) is None

