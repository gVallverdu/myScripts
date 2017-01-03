#!/usr/bin/env python

"""
Resize all images according to a given width.
Look at ImageMagick if available on your OS.
"""

__author__ = "Germain Salvato Vallverdu"

import os
import shutil
import argparse

try:
    from PIL import Image
except ImportError:
    print("The module PIL is needed !! sorry")
    raise ImportError("Need module PIL")


def list_img_files(folder="./", ext=(".jpg", ".jpeg", ".png", ".gif")):
    """ Return a list of image files according to ext list """
    return [f for f in os.listdir(folder) if any(f.lower().endswith(e) for e in ext)]


def get_options():
    """ get options from command lines """

    parser = argparse.ArgumentParser(description='Read options from the command line')
    parser.add_argument("-t", "--to",
                        help="final directory where resized images will be stored",
                        metavar="PATH",
                        default="./convert")
    parser.add_argument("-f", "--from",
                        help="working directory from where images to be resized have to be read",
                        metavar="PATH",
                        default="./")
    parser.add_argument("width",
                        help="Desired width of resized images",
                        metavar="WIDTH",
                        default=500,
                        type=int)

    return parser.parse_args()

if __name__ == "__main__":

    # manage options
    args = vars(get_options())  # convert to dict
    final_width = args["width"]
    final_dir = args["to"]
    workdir = args["from"]

    if os.path.exists(final_dir):
        print("Directory '{}' exists. Remove or rename it first !".format(final_dir))
        answer = input("Delete it ? (y/n) ")
        n = 0
        while answer.strip() != "n" and answer.strip() != "y" and n != 4:
            print("Please say 'y' for Yes, or 'n' for No")
            answer = input("Delete it ? (y/n) ")
            n += 1

        if answer == "y":
            shutil.rmtree(final_dir)
            print("Remove {}".format(final_dir))
        else:
            print("Nothing, was done")
            print("Remove or rename your output dir first")
            exit(0)

    os.mkdir(final_dir)

    for fimg in list_img_files(workdir):
        print(fimg)

        # open image file
        img = Image.open(fimg)
        w, h = img.size

        # compute resize factor according to the width
        rfactor = final_width / w

        # new width, height and image
        nw, nh = int(rfactor * w), int(rfactor * h)
        rimg = img.resize((nw, nh))

        # save image to new dir
        rimg.save(os.path.join(final_dir, fimg))
