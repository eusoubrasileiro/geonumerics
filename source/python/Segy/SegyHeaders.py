#
# Full structure of segy in python
# EBCDIC + Binary Header + N*(Trace Header)
# N number o traces
# Minimum size of a segy 3200 + 400 + 240

from ctypes import Structure, c_int, c_float, c_ushort, c_short, c_char

class SegyEBCDIC(Structure):
    """
    Ebcdic or text header 3200 bytes
    """
    _fields_ = [("ebcdic", 3200*c_char)]

class SegyBH(Structure):
    """ 
    Segy binary header 400 bytes 
    """
    _fields_ = [ 
        ("jobid", c_short),
        ("lino", c_short),
        ("reno", c_short),
        ("ntrpr", c_short),
        ("nart", c_short),
        ("hdt", c_ushort),
        ("dto", c_ushort),
        ("hns", c_ushort),
        ("nso", c_ushort),
        ("format", c_short),
        ("fold", c_short),
        ("tsort", c_short),
        ("vscode", c_short),
        ("hsfs", c_short),
        ("hsfe", c_short),
        ("hslen", c_short),
        ("hstyp", c_short),
        ("schn", c_short),
        ("hstas", c_short),
        ("hstae", c_short),
        ("htatyp", c_short),
        ("hcorr", c_short),
        ("bgrcv", c_short),
        ("rcvm", c_short),
        ("mfeet", c_short),
        ("polyt", c_short),
        ("vpol", c_short),
        ("hunass", 170*c_short)
    ]

class SegyTH(Structure):
    """
    Segy trace  header 240 bytes 
    """
    _fields_ = [    
        ("tracl", c_int),
        ("tracr", c_int),
        ("fldr", c_int),
        ("tracf", c_int),
        ("ep", c_int),
        ("cdp", c_int),
        ("cdpt", c_int),
        ("trid", c_short),
        ("nvs", c_short),
        ("nhs", c_short),
        ("duse", c_short),
        ("offset", c_int),
        ("gelev", c_int),
        ("selev", c_int),
        ("sdepth", c_int),
        ("gdel", c_int),
        ("sdel", c_int),
        ("swdep", c_int),
        ("gwdep", c_int),
        ("scalel", c_short),
        ("scalco", c_short),
        ("sx", c_int),
        ("sy", c_int),
        ("gx", c_int),
        ("gy", c_int),
        ("counit", c_short),
        ("wevel", c_short),
        ("swevel", c_short),
        ("sut", c_short),
        ("gut", c_short),
        ("sstat", c_short),
        ("gstat", c_short),
        ("tstat", c_short),
        ("laga", c_short),
        ("lagb", c_short),
        ("delrt", c_short),
        ("muts", c_short),
        ("mute", c_short),
        ("ns", c_ushort),
        ("dt", c_ushort),
        ("gain", c_short),
        ("igc", c_short),
        ("igi", c_short),
        ("corr", c_short),
        ("sfs", c_short),
        ("sfe", c_short),
        ("slen", c_short),
        ("styp", c_short),
        ("stas", c_short),
        ("stae", c_short),
        ("tatyp", c_short),
        ("afilf", c_short),
        ("afils", c_short),
        ("nofilf", c_short),
        ("nofils", c_short),
        ("lcf", c_short),
        ("hcf", c_short),
        ("lcs", c_short),
        ("hcs", c_short),
        ("year", c_short),
        ("day", c_short),
        ("hour", c_short),
        ("minute", c_short),
        ("sec", c_short),
        ("timbas", c_short),
        ("trwf", c_short),
        ("grnors", c_short),
        ("grnofr", c_short),
        ("grnlof", c_short),
        ("gaps", c_short),
        ("otrav", c_short),
        ("d1", c_float),
        ("f1", c_float),
        ("d2", c_float),
        ("f2", c_float),
        ("ungpow", c_float),
        ("unscale", c_float),
        ("ntr", c_short),
        ("mark", c_short),
        ("shortpad", c_short),
        ("unass", 14*c_short),
    ]
