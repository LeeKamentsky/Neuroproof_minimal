/* loadingplan.h - read a microns pipeline loading plan
 *
 * A loading plan is a .json file that gives the locations of the
 * TIFF stacks that must be loaded and composited to make a volume.
 * The format of a loading plan is described here:
 *
 * https://github.com/microns-ariadne/pipeline_engine/blob/master/README.md#load-plan-file-format
 *
 */
#ifndef _LOADINGPLAN
#define _LOADINGPLAN

#include <stddef.h>
#include <tiffio.h>
#include <json/value.h>
#include <json/reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cilk/cilk.h>
#include "ErrMsg.h"

/*
 * read_loading_plan - read the volume described by a loading plan into
 *                     a 3-D array
 * path -   path to the loading plan file
 * result - a pointer which will receive a 3-d array. The order is X, then Y
 *          then Z.
 * width -  the width (x extent) of the volume is stored here
 * height - the height (y extent) of the volume is stored here
 * depth -  the depth (z extent) of the volume is stored here
 */
template <typename T> void read_loading_plan(
    std::string path,
    T *&result,
    size_t &width,
    size_t &height,
    size_t &depth)
{
    Json::Reader reader;
    Json::Value d;
    
    std::cout << "Reading load plan, " << path << std::endl;
    ifstream fin(path);
    if (! fin) {
        throw NeuroProof::ErrMsg("Error: input file, \"" + path + "\" cannot be opened.");
    }
    if (! reader.parse(fin, d)) {
        throw NeuroProof::ErrMsg("Cannot parse \"" + path + "\" as json.");
    }
    fin.close();
    Json::Value dimensions = d["dimensions"];
    depth = (size_t)(dimensions[0].asUInt());
    height = (size_t)(dimensions[1].asUInt());
    width = (size_t)(dimensions[2].asUInt());
    Json::UInt x0 = d["x"].asUInt();
    Json::UInt y0 = d["y"].asUInt();
    Json::UInt z0 = d["z"].asUInt();
    Json::UInt x1 = x0 + width;
    Json::UInt y1 = y0 + height;
    Json::UInt z1 = z0 + depth;
    std::cout << "Total volume: x=" << x0 << " y=" << y0 << " z=" << z0;
    std::cout << " width=" << width << " height=" << height << " depth=" << depth << std::endl;
    result = new T[width * height * depth];
    TIFFErrorHandler old_handler = TIFFSetWarningHandler(NULL);

    Json::Value blocks = d["blocks"];
    cilk_for (int i=0; i < blocks.size(); i++) {
        Json::Value subvolume = blocks[i][1];
        Json::UInt svx0 = subvolume["x"].asUInt();
        Json::UInt svy0 = subvolume["y"].asUInt();
        Json::UInt svz0 = subvolume["z"].asUInt();
        Json::UInt svx1 = svx0 + subvolume["width"].asUInt();
        Json::UInt svy1 = svy0 + subvolume["height"].asUInt();
        Json::UInt svz1 = svz0 + subvolume["depth"].asUInt();
        std::string location = blocks[i][0].asString();
        std::cout << "Reading block from " << location << std::endl;
        std::cout << "   x=" << svx0 << ":" << svx1 << std::endl;
        std::cout << "   y=" << svy0 << ":" << svy1 << std::endl;
        std::cout << "   z=" << svz0 << ":" << svz1 << std::endl;
        if ((svx0 < x0) || (svx1 > x1) || (svy0 < y0) || (svy1 > y1) || (svz0 < z0) || (svz1 > z1)) {
            throw NeuroProof::ErrMsg("Sub-block out of bounds");
        }
        TIFF *tiff = TIFFOpen(location.c_str(), "r");
        if (! tiff) {
            throw NeuroProof::ErrMsg("Failed to read \"" + location + "\".");
        }
       
        /*
         * Special case: Python tifffile will take a volume that has 
         * a width or height of 1 and turn it into a single XZ or YZ plane
         *
         */
        if ((svx1 - svx0 == 1) || (svy1 - svy0 == 1)) {
            uint32 tiffWidth;
            uint32 tiffHeight;
            uint16 tiffBitsPerSample;
            uint16 tiffSamplesPerPixel;
            TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &tiffWidth);
            TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &tiffHeight);
            if (tiffWidth != (svx1 - svx0) * (svy1 - svy0)) {
                TIFFClose(tiff);
                std::cerr << "TIFF size not equal to expected" << std::endl;
                throw NeuroProof::ErrMsg("TIFF height not equal to expected");
            }
            if (tiffHeight != svz1 - svz0) {
                TIFFClose(tiff);
                std::cerr << "TIFF size not equal to expected" << std::endl;
                throw NeuroProof::ErrMsg("TIFF height not equal to expected");
            }
            TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &tiffBitsPerSample);
            if (tiffBitsPerSample / 8 != sizeof(T)) {
                TIFFClose(tiff);
                std::cerr << "# of bits / pixel does not match" << std::endl;
                throw NeuroProof::ErrMsg("# of bits / pixel does not match");
            }
            TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &tiffSamplesPerPixel);
            if (tiffSamplesPerPixel != 1) {
                TIFFClose(tiff);
                std::cerr << "More than one sample / pixel" << std::endl;
                throw NeuroProof::ErrMsg("More than one sample / pixel");
            }
            /*
             * Have to read individual voxels and advance in Y
             */
            if (svx1 - svx0 == 1) {
                std::cout << "Reading single plane as Y scans" << std::endl;
                std::unique_ptr<tdata_t> buf{ new tdata_t [tiffWidth] };
                for (unsigned int tiffZ=svz0; tiffZ < svz1; tiffZ++) {
                    TIFFReadScanline(tiff, &*buf, tiffZ-svz0);
                    for (unsigned int tiffY=svy0; tiffY < svy1; tiffY++) {
                        size_t offset = ((tiffZ - z0) * height + 
                                         (tiffY-y0)) * width +
                                         svx0 - x0;
                        result[offset] = ((T *)&*buf)[tiffY - svy0];
                    }
                }   
            } else {
                /* Read scan lines, but advance in Z */
                std::cout << "Reading single plane as X scans" << std::endl;
                for (unsigned int tiffZ=svz0; tiffZ < svz1; tiffZ++) {
                    size_t offset = ((tiffZ - z0) * height + (svy0-y0)) * width +
                                 svx0 - x0;
                    TIFFReadScanline(tiff, (tdata_t)(result+offset), tiffZ-svz0);
                }
            }
        } else {
            
            for (unsigned int tiffZ=svz0; tiffZ < svz1; tiffZ++) {
                uint32 tiffWidth;
                uint32 tiffHeight;
                uint16 tiffBitsPerSample;
                uint16 tiffSamplesPerPixel;
                TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &tiffWidth);
                if (tiffWidth != svx1 - svx0) {
                    TIFFClose(tiff);
                    std::cerr << "TIFF width not equal to expected" << std::endl;
                    throw NeuroProof::ErrMsg("TIFF width not equal to expected");
                }
                TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &tiffHeight);
                if (tiffHeight != svy1 - svy0) {
                    TIFFClose(tiff);
                    std::cerr << "TIFF height not equal to expected" << std::endl;
                    throw NeuroProof::ErrMsg("TIFF height not equal to expected");
                }
                TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &tiffBitsPerSample);
                if (tiffBitsPerSample / 8 != sizeof(T)) {
                    TIFFClose(tiff);
                    std::cerr << "# of bits / pixel does not match" << std::endl;
                    throw NeuroProof::ErrMsg("# of bits / pixel does not match");
                }
                TIFFGetField(tiff, TIFFTAG_SAMPLESPERPIXEL, &tiffSamplesPerPixel);
                if (tiffSamplesPerPixel != 1) {
                    TIFFClose(tiff);
                    std::cerr << "More than one sample / pixel" << std::endl;
                    throw NeuroProof::ErrMsg("More than one sample / pixel");
                }
               
                for (unsigned int tiffY=svy0; tiffY < svy1; tiffY++) {
                    size_t offset = ((tiffZ - z0) * height + (tiffY-y0)) * width +
                                 svx0 - x0;
                    TIFFReadScanline(tiff, (tdata_t)(result+offset), tiffY-svy0);
                }
                if (! TIFFReadDirectory(tiff)) {
                    if (tiffZ != svz1 - 1) {
                        std::cerr << "TIFF stack " << location << " contained too few slices" << std::endl;
                        break;
                    }
                }
            }
       }
       TIFFClose(tiff);
    }
    TIFFSetWarningHandler(old_handler);
    std::cout << "Finished reading volume." << std::endl;
}
#endif