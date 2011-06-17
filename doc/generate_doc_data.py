import os

outgd = " -o generated_data"
pipegd = " > generated_data/%s"

def probe_layouts():
    for probe in ["mar32.probe","mar16.probe","tetrode.probe"]:
        os.system("../plot_probe.py ../probes/%s"%probe + outgd)

def cmd_line_help():
    for script in ["cluster_from_raw_data.py","generalize_from_raw_data.py","check_crosstalk.py","plot_probe.py"]:
        os.system("../%s -h"%script + pipegd%(script[:-2] + "help.txt"))

def ch_colormaps():
    for datfile,probefile in [("marn2dec13/marn2dec13.002-5/marn2dec13.002-5.dat","aux/mar32_before.probe"),
                              ("marn2dec13/marn2dec13.002-5/marn2dec13.002-5.dat","aux/mar32_after.probe")]:
        os.system("../check_crosstalk.py /home/joschu/Data/mariano/%s -p %s"%(datfile,probefile)+outgd)

def features():
    os.system("../plot_features.py -s" + outgd)        

if __name__ == "__main__":
    probe_layouts()
    cmd_line_help()
    ch_colormaps()
    features()




