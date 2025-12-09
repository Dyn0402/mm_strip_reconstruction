//
// Created by dn277127 on 2025-12-02.
//

#ifndef dreamdataline_h
#define dreamdataline_h

inline bool is_final_trailer( uint16_t data )   { return (((data) & 0x7000)>>12)==7; }  // X111
inline bool is_data_trailer( uint16_t data )   { return (((data) & 0x6000)>>13)==2; }  // X10X
inline bool is_first_line( uint16_t data )   { return (((data) & 0x7000)>>12)==3; }  // X011
//inline bool is_data( uint16_t data )   { return (((data) & 0x7000)>>12)==0; }  // X000
inline bool is_data( uint16_t data )   { return (((data) & 0x6000)>>13)==0; }  // X00X
inline bool is_channel_ID( uint16_t data )   { return (((data) & 0x7000)>>12)==1; }  // X001
inline bool is_Feu_header( uint16_t data )   { return (((data) & 0x7000)>>12)==6; }  // X110
inline bool is_data_header( uint16_t data )   { return (((data) & 0x6000)>>13)==1; }  // X01X
inline bool get_zs_mode( uint16_t data )   { return (((data) & 0x400)>>10); }
inline uint16_t get_Feu_ID( uint16_t data )   { return (((data) & 0xFF)); };
inline uint16_t get_Event_ID( uint16_t data )   { return (((data) & 0xFFF)); }
inline uint16_t get_sample_ID( uint16_t data )   { return (((data) & 0xFF8)>>3); }
inline uint16_t get_channel_ID( uint16_t data )   { return (((data) & 0x3F)); }
inline uint16_t get_dream_ID( uint16_t data )   { return (((data) & 0xE00)>>9); }
inline uint16_t get_data( uint16_t data )   { return (((data) & 0xFFF)); }

//inline uint16_t get_dream_ID( uint16_t data )   { return (((data) & 0xE00)>>9); }
inline uint16_t get_dream_ID_ZS( uint16_t data )   { return (((data) & 0x1C0)>>6); }

inline uint16_t get_EoE( uint16_t data ){ return (((data) & 0x800)>>11); }
inline uint16_t get_timestamp( uint16_t data )   { return (((data) & 0xFFF)); }
inline uint16_t get_fine_timestamp( uint16_t data )   { return (((data) & 0x7)); }
#endif
