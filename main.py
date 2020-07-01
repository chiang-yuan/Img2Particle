import argparse

version = '1.0'

parser = argparse.ArgumentParser(prog='Img2Particle',
                                 description='Convert images into particle models.')

parser.add_argument('infile', metavar='images', type=str,
                    help='absolute or relative address of image files')
parser.add_argument('-o', dest='outfile', metavar='data file', type=argparse.FileType('w'))

parser.add_argument('--version', action='version', version='%(prog)s {:s}'.format(version))




args = parser.parse_args()
print(args)
