using UnityEngine;
using UnityEngine.Rendering;
using System;
using System.IO;

namespace BrunetonsImprovedAtmosphere
{
    static public class CBWrite
    {
        private static string[,] writeNames2D = new string[,]
        {
            {"write2DC1", "_Des2DC1", "_Buffer2DC1"},
            {"write2DC2", "_Des2DC2", "_Buffer2DC2"},
            {"write2DC3", "_Des2DC3", "_Buffer2DC3"},
            {"write2DC4", "_Des2DC4", "_Buffer2DC4"}
        };

        private static string[,] writeNames3D = new string[,]
        {
            {"write3DC1", "_Des3DC1", "_Buffer3DC1"},
            {"write3DC2", "_Des3DC2", "_Buffer3DC2"},
            {"write3DC3", "_Des3DC3", "_Buffer3DC3"},
            {"write3DC4", "_Des3DC4", "_Buffer3DC4"}
        };

        public static void IntoRenderTexture(RenderTexture tex, int channels, ComputeBuffer buffer, ComputeShader write)
        {
            Check(tex, channels, buffer, write);

            int kernel = -1;
            int depth = 1;

            if (tex.dimension == TextureDimension.Tex3D)
            {
                depth = tex.volumeDepth;
                kernel = write.FindKernel(writeNames3D[channels - 1, 0]);
                write.SetTexture(kernel, writeNames3D[channels - 1, 1], tex);
                write.SetBuffer(kernel, writeNames3D[channels - 1, 2], buffer);
            }
            else
            {
                kernel = write.FindKernel(writeNames2D[channels - 1, 0]);
                write.SetTexture(kernel, writeNames2D[channels - 1, 1], tex);
                write.SetBuffer(kernel, writeNames2D[channels - 1, 2], buffer);
            }

            if (kernel == -1)
                throw new ArgumentException("Could not find kernel " + writeNames2D[channels - 1, 0]);

            int width = tex.width;
            int height = tex.height;

            //set the compute shader uniforms
            write.SetInt("_Width", width);
            write.SetInt("_Height", height);
            write.SetInt("_Depth", depth);
            //run the  compute shader. Runs in threads of 8 so non divisable by 8 numbers will need
            //some extra threadBlocks. This will result in some unneeded threads running 
            int padX = (width % 8 == 0) ? 0 : 1;
            int padY = (height % 8 == 0) ? 0 : 1;
            int padZ = (depth % 8 == 0) ? 0 : 1;

            write.Dispatch(kernel, Mathf.Max(1, width / 8 + padX), Mathf.Max(1, height / 8 + padY), Mathf.Max(1, depth / 8 + padZ));
        }

        public static void IntoRenderTexture(RenderTexture tex, int channels, string Path, ComputeBuffer buffer, ComputeShader write)
        {
            Check(tex, channels, buffer, write);

            int kernel = -1;
            int depth = 1;

            if (tex.dimension == TextureDimension.Tex3D)
            {
                depth = tex.volumeDepth;
                kernel = write.FindKernel(writeNames3D[channels - 1, 0]);
                write.SetTexture(kernel, writeNames3D[channels - 1, 1], tex);
                write.SetBuffer(kernel, writeNames3D[channels - 1, 2], buffer);
            }
            else
            {
                kernel = write.FindKernel(writeNames2D[channels - 1, 0]);
                write.SetTexture(kernel, writeNames2D[channels - 1, 1], tex);
                write.SetBuffer(kernel, writeNames2D[channels - 1, 2], buffer);
            }

            if (kernel == -1)
                throw new ArgumentException("Could not find kernel " + writeNames2D[channels - 1, 0]);

            int width = tex.width;
            int height = tex.height;
            int size = width * height * depth * channels;

            float[] map = new float[size];
            LoadRawFile(Path, map, size);

            buffer.SetData(map);

            //set the compute shader uniforms
            write.SetInt("_Width", width);
            write.SetInt("_Height", height);
            write.SetInt("_Depth", depth);
            //run the  compute shader. Runs in threads of 8 so non divisable by 8 numbers will need
            //some extra threadBlocks. This will result in some unneeded threads running 
            int padX = (width % 8 == 0) ? 0 : 1;
            int padY = (height % 8 == 0) ? 0 : 1;
            int padZ = (depth % 8 == 0) ? 0 : 1;

            write.Dispatch(kernel, Mathf.Max(1, width / 8 + padX), Mathf.Max(1, height / 8 + padY), Mathf.Max(1, depth / 8 + padZ));
        }

        private static void LoadRawFile(string Path, float[] map, int size)
        {
            FileInfo fi = new FileInfo(Path);

            if (fi == null)
                throw new ArgumentException("Raw file not found (" + Path + ")");

            FileStream fs = fi.OpenRead();

            byte[] data = new byte[fi.Length];
            fs.Read(data, 0, (int)fi.Length);
            fs.Close();

            //divide by 4 as there are 4 bytes in a 32 bit float
            if (size > fi.Length / 4)
                throw new ArgumentException("Raw file is not the required size (" + Path + ")");

            for (int x = 0, i = 0; x < size; x++, i += 4)
            {
                //Convert 4 bytes to 1 32 bit float
                map[x] = System.BitConverter.ToSingle(data, i);
            };

        }

        private static void Check(RenderTexture tex, int channels, ComputeBuffer buffer, ComputeShader writeData)
        {
            if (tex == null)
                throw new ArgumentException("RenderTexture is null");

            if (buffer == null)
                throw new ArgumentException("Buffer is null");

            if (writeData == null)
                throw new ArgumentException("Computer shader is null");

            if (channels < 1 || channels > 4)
                throw new ArgumentException("Channels must be 1, 2, 3, or 4");

            if (!tex.enableRandomWrite)
                throw new ArgumentException("You must enable random write on render texture");

            if (!tex.IsCreated())
                throw new ArgumentException("Tex has not been created (Call Create() on tex)");
        }
    }
}




















