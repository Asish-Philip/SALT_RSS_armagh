from astropy.io import fits
import sys

arg = sys.argv
date = arg[1]#input('Please enter date: ')
file = arg[2]#input('enter file no: ')
filename = 'mbxgpP' + date + file.zfill(4) + '.fits'

hdul = fits.open('/Users/asishphilipmonai/SALT/raw_data/'+date+'/product/'+filename)
# name = './'+ date +'/r'+file+'.fits'
# hdul = fits.open(name)

hdr = hdul[0].header

print('+---------- FITS INFO ----------+')
hdul.info()
print('+------------ HEADER -----------+')
print(repr(hdr))
#

#######################################################
##               To edit fits header
#######################################################
# data, header = fits.getdata('/Users/asishphilipmonai/SALT/raw_data/'+date+'/product/'+filename, header=True)
#
# hdul = fits.open('/Users/asishphilipmonai/SALT/raw_data/'+date+'/product/'+filename)
# hdul.info()
# print(repr(hdul[0].header))
#
# hdul[0].header['OBJECT'] = 'ARC               '
# hdul[0].header['OBSTYPE'] = 'ARC               '
# hdul[0].header['CCDTYPE'] = 'ARC               '
#
# # hdul[0].header['OBSTYPE'] = 'OBJECT            '
# # hdul[0].header['CCDTYPE'] = 'OBJECT            '
# print('#########################')
# print(repr(hdul[0].header))
# hdul.writeto(filename, overwrite=True)
