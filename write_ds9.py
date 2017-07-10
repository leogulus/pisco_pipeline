from astropy import units as u
import os
import sys

#------
if __name__ == "__main__":
    print 'Number of arguments:', len(sys.argv), 'arguments.'
    print 'Argument List:', str(sys.argv)

    

def write_ds9_region_sky(table, allcol, color, fname):
    """
    write_ds9_region_sky: write out a ds9 region file
    INPUT:
    - table: the table that we are interested in plotting
    - allcol: the column name for the celestial coordinates
    - color: the color for the region in ds9 (usually 'green', 'red', 'blue')
    - fname: the output filename (should be 'xxx.reg')
    OUTPUT:
    - create a fname reg file for ds9 to see on the celestial coordinates
    """
    dir = "region/"
    if not os.path.exists(dir):
        os.makedirs(dir)
    text_file = open(os.path.join(dir, fname), "w")
    text_file.write("# Region file format: CIAO version 1.0\n")
    for i in range(0, len(table)):
        ra = table[allcol][i].ra.to_string(unit=u.hr, sep=':')
        dec = table[allcol][i].dec.to_string(unit=u.degree, sep=':')
        text_file.write("j2000; circle(%s,%s,0.1') # color=%s\n" %
                        (ra, dec, color))
    text_file.close()
