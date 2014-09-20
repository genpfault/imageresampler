// resampler.h, Separable filtering image rescaler v2.21, Rich Geldreich - richgel99@gmail.com
// See unlicense.org text at the bottom of this file.
#ifndef RESAMPLER_H
#define RESAMPLER_H

#include <vector>
#include <memory>
#include <map>

#define RESAMPLER_DEFAULT_FILTER "lanczos4"

// float or double
typedef float Resample_Real;

class Resampler
{
public:
    struct Contrib
    {
        Resample_Real weight;
        unsigned short pixel;
    };

    struct Contrib_List
    {
        unsigned short n;
        Contrib* p;
    };

    enum Boundary_Op
    {
        BOUNDARY_WRAP = 0,
        BOUNDARY_REFLECT = 1,
        BOUNDARY_CLAMP = 2
    };

    // Filter accessors.
    static unsigned int GetFilterCount();
    static const char* GetFilterName( unsigned int filter_num );

    class ContribLists
    {
    public:
        // src_w/src_h - Input dimensions
        // dst_w/dst_h - Output dimensions
        // boundary_op - How to sample pixels near the image boundaries
        // src_x_ofs/src_y_ofs - Offset input image by specified amount (fractional values okay)
        ContribLists
            (
            unsigned int src_w,
            unsigned int src_h,
            unsigned int dst_w,
            unsigned int dst_h,
            const char* Pfilter_name = RESAMPLER_DEFAULT_FILTER,
            Resample_Real filter_x_scale = 1.0f,
            Resample_Real filter_y_scale = 1.0f,
            Boundary_Op boundary_op = BOUNDARY_CLAMP,
            Resample_Real src_x_ofs = 0.0f,
            Resample_Real src_y_ofs = 0.0f
            );

        unsigned int GetSrcW() const { return m_Pclistc_x.srcDim; }
        unsigned int GetSrcH() const { return m_Pclistc_y.srcDim; }
        unsigned int GetDstW() const { return m_Pclistc_x.clists.size(); }
        unsigned int GetDstH() const { return m_Pclistc_y.clists.size(); }

        const Contrib_List* GetClistX() const { return &m_Pclistc_x.clists[ 0 ]; }
        const Contrib_List* GetClistY() const { return &m_Pclistc_y.clists[ 0 ]; }

    private:
        ContribLists();
        ContribLists( const ContribLists& );
        ContribLists& operator=( const ContribLists& );

        struct Contrib_List_Container
        {
            unsigned int srcDim;
            std::vector< Contrib > cpool;
            std::vector< Contrib_List > clists;
        };
        Contrib_List_Container m_Pclistc_x;
        Contrib_List_Container m_Pclistc_y;
    };

    typedef Resample_Real Sample;

    // Pclist_x/Pclist_y - Optional pointers to contributor lists from another instance of a Resampler
    // sample_low/sample_high - Clamp output samples to specified range, or disable clamping if sample_low >= sample_high
    Resampler
        (
        const ContribLists& contribLists,
        Resample_Real sample_low = 0.0f,
        Resample_Real sample_high = 0.0f
		);

    bool StartResample
        (
        unsigned int dst_subrect_x = 0,
        unsigned int dst_subrect_y = 0,
        unsigned int dst_subrect_w = 0,
        unsigned int dst_subrect_h = 0
        );

    void PutLine( const Sample* Psrc );

    // false if no scanlines are currently available (give the resampler more scanlines!)
    bool GetLine( Sample* Pdst );

private:
    Resampler();
    Resampler( const Resampler& );
    Resampler& operator=( const Resampler& );

    void resample_x( Sample* Pdst, const Sample* Psrc );
    void resample_y( Sample* Pdst );

    unsigned int m_intermediate_x;

    const unsigned int m_resample_src_w;
    const unsigned int m_resample_src_h;
    const unsigned int m_resample_dst_w;
    const unsigned int m_resample_dst_h;
   
    unsigned int m_dst_subrect_beg_x;
    unsigned int m_dst_subrect_end_x;
    unsigned int m_dst_subrect_beg_y;
    unsigned int m_dst_subrect_end_y;

    std::vector< Sample > m_Ptmp_buf;

    const Contrib_List* m_Pclist_x;
    const Contrib_List* m_Pclist_y;

    bool m_delay_x_resample;

    std::vector< int > m_Psrc_y_count_reference;
    std::vector< int > m_Psrc_y_count;

    std::vector< bool > m_Psrc_y_flag;
    std::map< int, std::vector< Sample > > m_Pscan_buf;

    unsigned int m_cur_src_y;
    unsigned int m_cur_dst_y;

    Resample_Real m_lo;
    Resample_Real m_hi;
};

#endif // RESAMPLER_H

// This is free and unencumbered software released into the public domain.
//
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
//
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//
// For more information, please refer to <http://unlicense.org/>
