/* storage_plan.h - write a volume according to a storage plan
 *
 * This implements writing a volume to a series of TIFF files according
 * to an ariadne_microns_pipeline storage plan (see link below for
 * the storage plan format):
 *
 * https://github.com/microns-ariadne/pipeline_engine/blob/master/README.md#storage-plan-file-format
 */

#ifndef _STORAGE_PLAN
#define _STORAGE_PLAN

#include <tiffio.h>
#include <string>
#include <json/value.h>
#include <json/reader.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cilk/cilk.h>
#include "ErrMsg.h"

using NeuroProof::ErrMsg;
/*
 * write_storage_plan - write a volume based on a storage plan
 *
 * path - path to the storage plan .json file
 * data - write the data from this volume. The data should be a raster in
 *        x, then y, then z of dimensions equal to width * height * depth
 * width - the width (x) of the data. This must correspond to the width
 *         of the volume in the storage plan
 * height - the height (y) of the data volume
 * depth - the depth (z) of the data volume
 */
template <typename T> void write_storage_plan(
    std::string path,
    T *data,
    size_t width,
    size_t height,
    size_t depth)
{
    Json::Reader reader;
    Json::Value d;

    std::cout << "Reading storage plan, " << path << std::endl;
    ifstream fin(path);
    if (! fin) {
        throw ErrMsg("Error: input file, \"" + path + "\" cannot be opened.");
    }
    if (! reader.parse(fin, d)) {
        throw ErrMsg("Cannot parse \"" + path + "\" as json.");
    }
    fin.close();
    Json::Value dimensions = d["dimensions"];
    Json::UInt jdepth = dimensions[0].asUInt();
    Json::UInt jheight = dimensions[1].asUInt();
    Json::UInt jwidth = dimensions[2].asUInt();
    if ((jwidth != width) || (jheight != height) || (jdepth != depth)) {
        throw ErrMsg("Volume dimensions differ from those in storage plan");
    }
    Json::UInt x0 = d["x"].asUInt();
    Json::UInt y0 = d["y"].asUInt();
    Json::UInt z0 = d["z"].asUInt();
    
    Json::Value blocks = d["blocks"];
    cilk_for (int i=0; i < blocks.size(); i++) {
        Json::Value subvolume = blocks[i][0];
        Json::UInt width = subvolume["width"].asUInt();
        Json::UInt height = subvolume["height"].asUInt();
        Json::UInt depth = subvolume["depth"].asUInt();
        Json::UInt svx0 = subvolume["x"].asUInt();
        Json::UInt svy0 = subvolume["y"].asUInt();
        Json::UInt svz0 = subvolume["z"].asUInt();
        Json::UInt svx1 = svx0 + width;
        Json::UInt svy1 = svy0 + height;
        Json::UInt svz1 = svz0 + depth;
        std::string location = blocks[i][1].asString();
        cout << "Writing " << location << endl;
        cout << "  x=" << svx0 << ":" << svx1;
        cout << "  y=" << svy0 << ":" << svy1;
        cout << "  z=" << svz0 << ":" << svz1;
        TIFF *tiff = TIFFOpen(location.c_str(), "w");
        if (! tiff) {
            throw ErrMsg("Unable to open " + location + " for writing.");
        }
        try {
            for (Json::UInt tiffZ=svz0; tiffZ < svz1; tiffZ++) {
                TIFFSetField(tiff, TIFFTAG_IMAGEWIDTH, width);
                TIFFSetField(tiff, TIFFTAG_IMAGELENGTH, height);
                TIFFSetField(tiff, TIFFTAG_BITSPERSAMPLE,
                             sizeof(T) * 8);
                TIFFSetField(tiff, TIFFTAG_SAMPLESPERPIXEL, 1);
                TIFFSetField(tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
                TIFFSetField(tiff, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
                TIFFSetField(tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
                TIFFSetField(tiff, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS);
                for (Json::UInt tiffY=svy0; tiffY < svy1; tiffY++) {
                    size_t offset = ((tiffZ - z0) * height + (tiffY-y0)) * width +
                                    svx0 - x0;
                    TIFFWriteScanline(
                        tiff, (void *)(data + offset), tiffY-svy0, 0);
                }
                TIFFWriteDirectory(tiff);
            }
        } catch (...) {
            TIFFClose(tiff);
            throw;
        }
        TIFFClose(tiff);
    }
    std::cout << "Finished writing volume" << std::endl;
}

#endif