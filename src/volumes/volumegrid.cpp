
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// volumes/volumegrid.cpp*
#include "stdafx.h"
#include "volumes/volumegrid.h"
#include "paramset.h"
#include <iostream>
#include <fstream>
#include <stdio.h>

using namespace std;

// VolumeGridDensity Method Definitions
float VolumeGridDensity::Density(const Point &Pobj) const {
    if (!extent.Inside(Pobj)) return 0;
    // Compute voxel coordinates and offsets for _Pobj_
    Vector vox = extent.Offset(Pobj);
    vox.x = vox.x * nx - .5f;
    vox.y = vox.y * ny - .5f;
    vox.z = vox.z * nz - .5f;
    int vx = Floor2Int(vox.x), vy = Floor2Int(vox.y), vz = Floor2Int(vox.z);
    float dx = vox.x - vx, dy = vox.y - vy, dz = vox.z - vz;

    // Trilinearly interpolate density values to compute local density
    float d00 = Lerp(dx, D(vx, vy, vz),     D(vx+1, vy, vz));
    float d10 = Lerp(dx, D(vx, vy+1, vz),   D(vx+1, vy+1, vz));
    float d01 = Lerp(dx, D(vx, vy, vz+1),   D(vx+1, vy, vz+1));
    float d11 = Lerp(dx, D(vx, vy+1, vz+1), D(vx+1, vy+1, vz+1));
    float d0 = Lerp(dy, d00, d10);
    float d1 = Lerp(dy, d01, d11);
    return Lerp(dz, d0, d1);
}


void LoadHeader(char *prefix,
                int &volumeWidth, int &volumeHeight, int &volumeDepth)
{
    char hdrFile[300];
    std::ifstream inputFileStream;

    // Adding the ".hdr" prefix to the dataset path
    sprintf(hdrFile, "%s.hdr", prefix);

    // Open the eader hdrFile to read the dataset dimensions
    inputFileStream.open(hdrFile, std::ios::in);

    // Checking the openning of the file
    if (inputFileStream.fail())
    {
        printf("Could not open the HDR file :%s", hdrFile);
        exit(0);
    }

    // Reading the dimensions
    inputFileStream >> volumeWidth;
    inputFileStream >> volumeHeight;
    inputFileStream >> volumeDepth;

    // Closing the ".hdr" file
    inputFileStream.close();
}

float* LoadVolume(char* prefix)
{
    char imgFile[100];
    ifstream inputFileStream;
    int iWidth;
    int iHeight;
    int iDepth;

    // Reading the header file
    LoadHeader(prefix, iWidth, iHeight, iDepth);

    // printf("Volume size of %d x %d x %d \n", iWidth, iHeight, iDepth);

    // Adding the ".img" prefix to the dataset path
    sprintf(imgFile, "%s.img", prefix);

    // Reading the volume image (luminance values)
    inputFileStream.open(imgFile, ios::in);
    if (inputFileStream.fail())
    {
        printf("Could not open %s \n", imgFile);
        return NULL;
    }

    // Allocating the luminance image
    char* dataChar = (char*) malloc (sizeof(char) * iWidth * iHeight * iDepth);
    float* data = (float*) malloc (sizeof(float) * iWidth * iHeight * iDepth);

    // Read the image byte by byte
    const int numVoxels = iWidth * iHeight * iDepth;
    inputFileStream.read(((char *)dataChar), numVoxels);

    // Closing the input volume stream
    inputFileStream.close();

    for (int i = 0; i < numVoxels; i++)
    {
        if (dataChar[i] > 100)
            data[i] = (float) 0;
        else
            data[i] = (float) dataChar[i];
    }

    free (dataChar);

    return data;
}

VolumeGridDensity *CreateGridVolumeRegion(const Transform &volume2world,
        const ParamSet &params) {

    // Initialize common volume region parameters
    Spectrum sigma_a = params.FindOneSpectrum("sigma_a", 0.);
    Spectrum sigma_s = params.FindOneSpectrum("sigma_s", 0.);
    float g = params.FindOneFloat("g", 0.);
    Spectrum Le = params.FindOneSpectrum("Le", 0.);
    Point p0 = params.FindOnePoint("p0", Point(0,0,0));
    Point p1 = params.FindOnePoint("p1", Point(1,1,1));

    int nx = params.FindOneInt("nx", 1);
    int ny = params.FindOneInt("ny", 1);
    int nz = params.FindOneInt("nz", 1);

    printf("NX=%d, NY=%d, NZ=%d \n", nx, ny, nz);

    int nitems = 10;
    const float *data = params.FindFloat("density", &nitems);
    if (!data) {
        Error("No \"density\" values provided for volume grid?");
        return NULL;
    }

    if (nitems != nx*ny*nz) {
        Error("VolumeGridDensity has %d density values but nx*ny*nz = %d",
            nitems, nx*ny*nz);
        return NULL;
    }

    // IF LOADING FROM A RAW FILE
    // char path[120] = "/PBR/Software/Scenes/RealStack/G_Volume/G_Volume";
    // const float* data = LoadVolume(path);

    return new VolumeGridDensity(sigma_a, sigma_s, g, Le, BBox(p0, p1),
        volume2world, nx, ny, nz, data);
}


