// resampler test, Rich Geldreich - richgel99@gmail.com
// See unlicense.org text at the bottom of resampler.h
// Example usage: resampler.exe input.tga output.tga width height
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <algorithm>

#include "resampler.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int main(int arg_c, char** arg_v)
{
   if (arg_c != 9 && arg_c != 5)
   {
      printf("Usage: input_image output_image.tga width height [subrect-x] [subrect-y] [subrect-width] [subrect-height]\n");
      return EXIT_FAILURE;
   }
   
   const char* pSrc_filename = arg_v[1];
   const char* pDst_filename = arg_v[2];
   const int dst_width = atoi(arg_v[3]);
   const int dst_height = atoi(arg_v[4]);

   const int subrect_x = arg_c == 9 ? atoi( arg_v[5] ) : 0;
   const int subrect_y = arg_c == 9 ? atoi( arg_v[6] ) : 0;
   const int subrect_w = arg_c == 9 ? atoi( arg_v[7] ) : dst_width;
   const int subrect_h = arg_c == 9 ? atoi( arg_v[8] ) : dst_height;
   
   if (std::min(dst_width, dst_height) < 1)
   {
      printf("Invalid output width/height!\n");
      return EXIT_FAILURE;
   }

   if ( subrect_w > dst_width )
   {
      printf("Invalid horizontal crop!\n");
      return EXIT_FAILURE;
   }

   if ( subrect_h > dst_height )
   {
      printf("Invalid vertical crop!\n");
      return EXIT_FAILURE;
   }
   
   printf("Loading image: %s\n", pSrc_filename);
   
   int src_width, src_height, n;
   unsigned char* pSrc_image = stbi_load(pSrc_filename, &src_width, &src_height, &n, 0);
   if (!pSrc_image)
   {
      printf("Failed loading image!\n");
      return EXIT_FAILURE;
   }
   printf("Resolution: %ux%u, Channels: %u\n", src_width, src_height, n);
   
   const int max_components = 4;   
   
   if (n > max_components)
   {
      printf("Image has too many components!\n");
      return EXIT_FAILURE;
   }
      
   // Partial gamma correction looks better on mips. Set to 1.0 to disable gamma correction. 
   const float source_gamma = 1.75f;
   
   // Filter scale - values < 1.0 cause aliasing, but create sharper looking mips.
   const float filter_scale = 1.0f;//.75f;
   
   const char* pFilter = "blackman";//RESAMPLER_DEFAULT_FILTER;
         
   float srgb_to_linear[256];
   for (int i = 0; i < 256; ++i)
      srgb_to_linear[i] = (float)pow(i * 1.0f/255.0f, source_gamma);

   const int linear_to_srgb_table_size = 4096;
   unsigned char linear_to_srgb[linear_to_srgb_table_size];
   
   const float inv_linear_to_srgb_table_size = 1.0f / linear_to_srgb_table_size;
   const float inv_source_gamma = 1.0f / source_gamma;

   for (int i = 0; i < linear_to_srgb_table_size; ++i)
   {
      int k = (int)(255.0f * pow(i * inv_linear_to_srgb_table_size, inv_source_gamma) + .5f);
      if (k < 0) k = 0; else if (k > 255) k = 255;
      linear_to_srgb[i] = (unsigned char)k;
   }
   
   Resampler::ContribLists* clists = new Resampler::ContribLists( src_width, src_height, dst_width, dst_height, pFilter, filter_scale, filter_scale, Resampler::BOUNDARY_CLAMP, 0.0f, 0.0f );

   Resampler* resamplers[max_components];
   std::vector<float> samples[max_components];

   // Now create a Resampler instance for each component to process. The first instance will create new contributor tables, which are shared by the resamplers 
   // used for the other components (a memory and slight cache efficiency optimization).
   for (int i = 0; i < n; i++)
   {
      resamplers[i] = new Resampler( *clists, 0.0f, 1.0f );
      resamplers[i]->StartResample( subrect_x, subrect_y, subrect_w, subrect_h );
      samples[i].resize(src_width);
   }      
      
   std::vector<unsigned char> dst_image( subrect_w * n * subrect_h );
   std::vector<float> pOutput_samples( subrect_w );

   const int src_pitch = src_width * n;
   const int dst_pitch = subrect_w * n;
   int dst_y = 0;
   
   printf("Resampling to %ux%u\n", dst_width, dst_height);
      
   for (int src_y = 0; src_y < src_height; src_y++)
   {
      const unsigned char* pSrc = &pSrc_image[src_y * src_pitch];
         
      for (int x = 0; x < src_width; x++)
      {
         for (int c = 0; c < n; c++)
         {
            if ((c == 3) || ((n == 2) && (c == 1)))
               samples[c][x] = *pSrc++ * (1.0f/255.0f);
            else
               samples[c][x] = srgb_to_linear[*pSrc++];        
         }
      }
      
      for (int c = 0; c < n; c++)         
      {
         resamplers[c]->PutLine(&samples[c][0]);
      }         
         
      for ( ; ; )
      {
         int comp_index;
         for (comp_index = 0; comp_index < n; comp_index++)
         {
            if( !resamplers[comp_index]->GetLine( &pOutput_samples[0] ) )
            {
               continue;
            }
            
            const bool alpha_channel = (comp_index == 3) || ((n == 2) && (comp_index == 1));
            assert(dst_y < dst_height);
            unsigned char* pDst = &dst_image[dst_y * dst_pitch + comp_index];
            
            for (int x = 0; x < subrect_w; x++)
            {
               if (alpha_channel)
               {
                  int c = (int)(255.0f * pOutput_samples[x] + .5f);
                  if (c < 0) c = 0; else if (c > 255) c = 255;
                  *pDst = (unsigned char)c;
               }
               else
               {
                  int j = (int)(linear_to_srgb_table_size * pOutput_samples[x] + .5f);
                  if (j < 0) j = 0; else if (j >= linear_to_srgb_table_size) j = linear_to_srgb_table_size - 1;
                  *pDst = linear_to_srgb[j];
               }
               
               pDst += n;
            }
         }     
         if (comp_index < n)
            break; 
         
         dst_y++;
      }
   }
   
   printf("Writing TGA file: %s\n", pDst_filename);
   
   if (!stbi_write_tga(pDst_filename, subrect_w, subrect_h, n, &dst_image[0]))
   {
      printf("Failed writing output image!\n");
      return EXIT_FAILURE;
   }
   
   stbi_image_free(pSrc_image);

   // Delete the resamplers.
   for (int i = 0; i < n; i++)
      delete resamplers[i];
   delete clists;
   
   return EXIT_SUCCESS;
}
