using UnityEngine;
using UnityEngine.Rendering;
using System;
using System.IO;

namespace BrunetonsImprovedAtmosphere
{
    static public class CBRead
    {

        private static string[,] readNames2D = new string[,]
        {
            {"read2DC1", "_Tex2D", "_Buffer2DC1"},
            {"read2DC2", "_Tex2D", "_Buffer2DC2"},
            {"read2DC3", "_Tex2D", "_Buffer2DC3"},
            {"read2DC4", "_Tex2D", "_Buffer2DC4"}
        };

        private static string[,] readNames3D = new string[,]
        {
            {"read3DC1", "_Tex3D", "_Buffer3DC1"},
            {"read3DC2", "_Tex3D", "_Buffer3DC2"},
            {"read3DC3", "_Tex3D", "_Buffer3DC3"},
            {"read3DC4", "_Tex3D", "_Buffer3DC4"}
        };

        public static void FromRenderTexture(RenderTexture tex, int channels, ComputeBuffer buffer, ComputeShader read)
        {
            Check(tex, channels, buffer, read);

            int kernel = -1;
            int depth = 1;

            if (tex.dimension == TextureDimension.Tex3D)
            {
                depth = tex.volumeDepth;
                kernel = read.FindKernel(readNames3D[channels - 1, 0]);
                read.SetTexture(kernel, readNames3D[channels - 1, 1], tex);
                read.SetBuffer(kernel, readNames3D[channels - 1, 2], buffer);
            }
            else
            {
                kernel = read.FindKernel(readNames2D[channels - 1, 0]);
                read.SetTexture(kernel, readNames2D[channels - 1, 1], tex);
                read.SetBuffer(kernel, readNames2D[channels - 1, 2], buffer);
            }

            if (kernel == -1)
                throw new ArgumentException("Could not find kernel " + readNames2D[channels - 1, 0]);

            int width = tex.width;
            int height = tex.height;

            read.SetInt("_Width", width);
            read.SetInt("_Height", height);
            read.SetInt("_Depth", depth);

            //run the  compute shader. Runs in threads of 8 so non divisable by 8 numbers will need
            //some extra threadBlocks. This will result in some unneeded threads running 
            int padX = (width % 8 == 0) ? 0 : 1;
            int padY = (height % 8 == 0) ? 0 : 1;
            int padZ = (depth % 8 == 0) ? 0 : 1;

            read.Dispatch(kernel, Mathf.Max(1, width / 8 + padX), Mathf.Max(1, height / 8 + padY), Mathf.Max(1, depth / 8 + padZ));

        }

        public static void SingleFromRenderTexture(RenderTexture tex, float x, float y, float z, ComputeBuffer buffer, ComputeShader read, bool useBilinear)
        {
            Check(tex, 0, buffer, read);

            int kernel = -1;
            int depth = 1;

            if (tex.dimension == TextureDimension.Tex3D)
            {
                if (useBilinear)
                    kernel = read.FindKernel("readSingleBilinear3D");
                else
                    kernel = read.FindKernel("readSingle3D");

                depth = tex.volumeDepth;
                read.SetTexture(kernel, "_Tex3D", tex);
                read.SetBuffer(kernel, "_BufferSingle3D", buffer);
            }
            else
            {
                if (useBilinear)
                    kernel = read.FindKernel("readSingleBilinear2D");
                else
                    kernel = read.FindKernel("readSingle2D");

                read.SetTexture(kernel, "_Tex2D", tex);
                read.SetBuffer(kernel, "_BufferSingle2D", buffer);
            }

            if (kernel == -1)
                throw new ArgumentException("Could not find kernel readSingle for " + tex.dimension);

            int width = tex.width;
            int height = tex.height;

            //used for point sampling
            read.SetInt("_IdxX", (int)x);
            read.SetInt("_IdxY", (int)y);
            read.SetInt("_IdxZ", (int)z);
            //used for bilinear sampling
            read.SetVector("_UV", new Vector4(x / (float)(width - 1), y / (float)(height - 1), z / (float)(depth - 1), 0.0f));

            read.Dispatch(kernel, 1, 1, 1);
        }

        private static void Check(RenderTexture tex, int channels, ComputeBuffer buffer, ComputeShader read)
        {
            if (tex == null)
                throw new ArgumentException("RenderTexture is null");

            if (buffer == null)
                throw new ArgumentException("Buffer is null");

            if (read == null)
                throw new ArgumentException("Computer shader is null");

            if (channels < 1 || channels > 4)
                throw new ArgumentException("Channels must be 1, 2, 3, or 4");

            if (!tex.IsCreated())
                throw new ArgumentException("Tex has not been created (Call Create() on tex)");
        }

    }
}




















