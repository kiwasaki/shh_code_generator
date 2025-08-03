#pragma once

#include<iostream>
#include<array>
#include<string>
#include<sstream>
#include<fstream>
#include<print>
#include<numbers>

#include"util_stable.hpp"

//Spherical Harmonicsを計算するコードを生成するクラス
// L: 球面調和関数の次数
template< size_t L > class shh_code_generator
{
public:

    //コンストラクタ
    explicit shh_code_generator( const std::string &prefix, const bool with_gradient = false, const bool with_gradient_hessian = false )
    {
        //ファイル名はデフォルトで sheval_00L.hpp
        const std::string filename = std::format( "{}_{:03}.hpp", prefix, L );
        std::ofstream fout( filename );
        if( fout.fail() ) {
            std::println( stderr, "cannot open file {}", filename );
            exit( - 1 );
        }

        if( with_gradient && not( with_gradient_hessian ) ) {
            fout << generate_gradient( prefix );
        } else if( with_gradient_hessian ) {
            fout << generate_hessian( prefix );
        } else {
            fout << generate( prefix );
        }
    }

private:
    //球面調和関数のid (l,m)->l*(l+1)+m
    [[ nodiscard ]] int id( const int l, const int m ) const { return l * ( l + 1 ) + m; }
    [[ nodiscard ]] int qlm_id( const int l, const int m ) const { return l * L + m; }

    //SHを計算するコードの生成
    std::string generate( const std::string &function_name = "sheval" ) const
    {
        using namespace util;

        std::ostringstream code;
        int l, m;
        //add comments
        code << "// x, y, z are normalized coordinates of the evaluation direction\n";
        code << "#pragma once\n";
        code << std::format( "template< size_t L > void {}( const float x, const float y, const float z, float *ylm );\n", function_name );
        code << std::format( "template<> void {}< {:3d} >( const float x, const float y, const float z, float *ylm )\n", function_name, L );
        code << "{\n";
        code << tab() << "float c0, s0, c1, s1, tmp0, tmp1, tmp2;\n";
        code << tab() << "float z2 = z * z;\n";
        code << tab() << "//zonal harmonics(m=0)\n";
        l = 0; m = 0;
        code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), constant( qmm( m ) ) );
        l = 1;
        code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), qm1m( m ) );
        if( L > 2 ) {
            l = 2;
            code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), qm2m( m ) );
        }
        for( l = 3; l < L; ++l ) {
            const auto coef0 = klm( l, m ) / klm( l - 1, m ) * ( 2.0 * l - 1.0 ) / ( l - m );
            const auto coef1 = klm( l, m ) / klm( l - 2, m ) * ( l + m - 1.0 ) / ( l - m );
            code << tab() << std::format( "ylm[ {:3d} ] = {} * z * ylm[ {:3d} ] - {} * ylm[ {:3d} ];\n", id( l, m ), constant( coef0 ), id( l - 1, m ), constant( coef1 ), id( l - 2, m ) );
        }
        code << tab() << std::format( "c0 = x; s0 = y;\n" );
        std::array< std::string, 2 > sc = { "c0", "c1" };
        std::array< std::string, 2 > ss = { "s0", "s1" };
        std::array< std::string, 3 > tmp = { "tmp0", "tmp1", "tmp2" };
        int idxsc = 0; //sine/cosine active index
        int idxp;
        //mについてのループ
        //1. y_m^m (y_m^{-m})を計算
        //2. y_{m+1}^mをy_m^mから計算
        //3. y_l^mをy_{l-1}^m, y_{l-2}^mから計算
        for( m = 1; m < L - 1; ++m ) {
            code << tab() << std::format( "//m = {:03d}\n", m );
            l = m;
            idxp = 0;
            code << tab() << std::format( "{} = {};\n", tmp[ idxp ], constant( qmm( m ) ) );
            code << tab() << std::format( "ylm[ {:3d} ] = {} * c0;\n", id( l, + m ), tmp[ idxp ], sc[ idxsc & 1 ] );
            code << tab() << std::format( "ylm[ {:3d} ] = {} * s0;\n", id( l, - m ), tmp[ idxp ], ss[ idxsc & 1 ] );

            if( m + 1 < L ) {
                l++;
                idxp++;
                code << tab() << std::format( "{} = {};\n", tmp[ idxp ], qm1m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * c0;\n", id( l, + m ), tmp[ idxp ] );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * s0;\n", id( l, - m ), tmp[ idxp ] );
            }

            if( m + 2 < L ) {
                l++;
                idxp++;
                code << tab() << std::format( "{} = {};\n", tmp[ idxp ], qm2m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * c0;\n", id( l, + m ), tmp[ idxp ] );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * s0;\n", id( l, - m ), tmp[ idxp ] );
            }

            for( l = m + 3; l < L; ++l ) {
                idxp++;
                code << tab() << std::format( "{} = {};\n", tmp[ idxp % 3 ], qlm( l, m, tmp[ ( idxp + 3 - 1 ) % 3 ], tmp[ ( idxp + 3 - 2 ) % 3 ] ) );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * c0;\n", id( l, + m ), tmp[ idxp % 3 ] );
                code << tab() << std::format( "ylm[ {:3d} ] = {} * s0;\n", id( l, - m ), tmp[ idxp % 3 ] );
            }
            code << tab() << "s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
        }
        //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        idxp = ( idxp + 1 ) % 3;
        code << tab() << std::format( "//m = {:03d}\n", m );
        code << tab() << std::format( "{} = {};\n", tmp[ idxp ], constant( qmm( m ) ) ); //assign( tmp[ idxp ], constant( qmm( m ) ) ) << ";\n";
        code << tab() << std::format( "ylm[ {:3d} ] = {} * c0;\n", id( l, + m ), tmp[ idxp ] );
        code << tab() << std::format( "ylm[ {:3d} ] = {} * s0;\n", id( l, - m ), tmp[ idxp ] );
        code << "}\n";
        return code.str();
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Spherical Harmonicsの空間的勾配を計算するコードの生成
    std::string generate_gradient( const std::string function_name = "shgrad" ) const
    {
        using namespace util;

        std::ostringstream code;
        int l, m;
        std::array< std::string, 4 > tmp = { "tmp0", "tmp1", "tmp2", "tmp3" };
        int idxp;

        //add comments
        code << "// x, y, z are the normalized coordinates of the evaluation direction\n";
        code << "#pragma once\n";
        code << std::format( "template< size_t L > void {}( const float x, const float y, const float z, float *ylm, std::array< glm::vec3, L * L > &glm );\n", function_name );
        code << std::format( "template<> void {}< {:3d} >( const float x, const float y, const float z, float *ylm, std::array< glm::vec3, {:3d} > &glm )\n", function_name, L, L * L );
        code << "{\n";
        code << tab() << "float c0, c1, s0, s1, tmp, tmp0, tmp1, tmp2, tmp3;\n";
        code << tab() << "float z2 = z * z;\n";
        code << tab() << std::format( "std::array< float, {:3d} > qlm {{}};\n", L * L );
        code << tab() << "//zonal harmonics(m=0)\n";
        l = 0; m = 0;
        code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), constant( qmm( m ) ) );
        code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
        l = 1;
        code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), qm1m( m ) );
        code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm1m( m ) );
        if( L > 2 ) {
            l = 2;
            code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), qm2m( m ) );
            code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm2m( m ) );
        }
        for( l = 3; l < L; ++l ) {
            code << tab() << std::format( "ylm[ {:3d} ] = {};\n", id( l, m ), qlm( l, m, ylm( l - 1, m ), ylm( l - 2, m ) ) ) ;
            code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qlm( l, m, std::format( "qlm[ {:3d} ]", qlm_id( l - 1, m ) ), std::format( "qlm[ {:3d} ]", qlm_id( l - 2, m ) ) ) );
        }
        code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0;\n";
        //mについてのループ
        //1. y_m^m (y_m^{-m})を計算
        //2. y_{m+1}^mをy_m^mから計算
        //3. y_l^mをy_{l-1}^m, y_{l-2}^mから計算
        for( m = 1; m < L - 1; ++m ) {
            code << tab() << std::format( "//m = {:03d}\n", m );
            l = m;
            idxp = 0;
            code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );

            if( m + 1 < L ) {
                l++;
                idxp++;
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm1m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }

            if( m + 2 < L ) {
                l++;
                idxp++;
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm2m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }

            for( l = m + 3; l < L; ++l ) {
                idxp++;
                //code << tab() << assign( tmp[ idxp % 3 ], qlm( l, m, tmp[ ( idxp + 3 - 1 ) % 3 ], tmp[ ( idxp + 3 - 2 ) % 3 ] ) ) << ";\n";
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qlm( l, m, std::format( "qlm[ {:3d} ]", qlm_id( l - 1, m ) ), std::format( "qlm[ {:3d} ]", qlm_id( l - 2, m ) ) ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }
            //update cosine and sine
            code << tab() << "s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
        }
        //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        idxp = ( idxp + 1 ) % 3;
        code << tab() << std::format( "//m = {:03d}\n", m );
        code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );

        //calculate gradient
        // m = 0 Zonal Harmonics
        code << tab() << "//calculate gradient\n";
        code << tab() << "const glm::vec3 zero = glm::vec3( 0.f, 0.f, 0.f );\n";
        l = m = 0;
        code << tab() << std::format( "glm[ {:3d} ] = zero;\n", id( l, m ) );
        l = 1;
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * ylm[ {0:3d} ];\n", id( l, m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * ylm[ {0:3d} ];\n", id( l, m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * ylm[ {0:3d} ] + {1};\n", id( l, m ), constant( klm( l, m ) ) );
        for( l = 2; l < L; ++l ) {
            code << tab() << std::format( "tmp0 = {:3d} * ylm[ {:3d} ]; tmp1 = {} * qlm[ {:3d} ];\n", l, id( l, m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ) );
            code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp0 - tmp1 );\n", id( l, m ) );
            code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp0 - tmp1 );\n", id( l, m ) );
            code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp0 + {} * qlm[ {:3d} ];\n", id( l, m ), constant( klm( l, m ) / klm( l - 1, m ) * l ), qlm_id( l - 1, m ) );
        }
        code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0;\n";

        for( m = 1; m < L - 1; ++m ) {
            code << tab() << std::format( "//m = {:03d}\n", m );
            for( l = m; l < L; ++l ) {
                code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c0; tmp1 = tmp * s0;\n", constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ) );
                code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp2 = tmp * c1; tmp3 = tmp * s1;\n", m, qlm_id( l, m ) );
                code << tab() << std::format( "tmp = {} * ylm[ {:3d} ];\n", l, id( l, + m ) );
                //grad x & y
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp - tmp0 ) + tmp2;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp - tmp0 ) - tmp3;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp;\n", id( l, + m ) );
                code << tab() << std::format( "tmp = {} * ylm[ {:3d} ];\n", l, id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp - tmp1 ) + tmp3;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp - tmp1 ) + tmp2;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp;\n", id( l, - m ) );
                code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; glm[ {:3d} ].z += tmp * c0; glm[ {:3d} ].z += tmp * s0;\n", constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), qlm_id( l - 1, m ), id( l, + m ), id( l, - m ) );
            }
            code << tab() << "s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
        }
        //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c1; tmp1 = tmp * s1; tmp2 = {} * ylm[ {:3d} ]; tmp3 = {} * ylm[ {:3d} ];\n", m, qlm_id( l, m ), l, id( l, + m ), l, id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * tmp2 + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * tmp2 - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * tmp2;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * tmp3 + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * tmp3 + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * tmp3;\n", id( l, - m ) );
        code << "}\n";
        return code.str();
    }

#pragma region generate_gradient_hessian_opt
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Spherical Harmonicsのhessianを計算するコードの生成 乗算数を減らすように工夫
    std::string generate_gradient_hessian_opt( const std::string function_name = "shhessian" ) const
    {
        using namespace util;

        std::ostringstream code;
        int l, m;
        //add comments
        code << "// x, y, z are normalized coordinates of the evaluation direction\n";
        code << "#pragma once\n";
        code << std::format( "template< size_t L > void {}( const float x, const float y, const float z, float *ylm, std::array< glm::vec3, L * L > &glm, std::array< glm::mat3, L * L > &hlm );\n", function_name );
        code << std::format( "template<> void {0}< {1:3d} >( const float x, const float y, const float z, float *ylm, std::array< glm::vec3, {2:3d} > &glm, std::array< glm::mat3, {2:3d} > &hlm )\n", function_name, L, L * L );
        code << "{\n";
        code << tab() << "float c0, c1, c2, s0, s1, s2, tmp, tmp0, tmp1, tmp2, tmp3;\n";
        code << tab() << "const float x2 = x * x;\n";
        code << tab() << "const float y2 = y * y;\n";
        code << tab() << "const float z2 = z * z;\n";
        code << tab() << "const float xy = x * y;\n";
        code << tab() << "const float xz = x * z;\n";
        code << tab() << "const float yz = y * z;\n";
        code << tab() << std::format( "std::array< float, {:3d} > qlm {{}};\n", L * L );
        code << tab() << "//zonal harmonics(m=0)\n";
        l = 0; m = 0;
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] = {};\n", id( l, m ), qlm_id( l, m ), constant( qmm( m ) ) );
        l = 1;
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] = {};\n", id( l, m ), qlm_id( l, m ), qm1m( m ) );
        if( L > 2 ) {
            l = 2;
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] = {};\n", id( l, m ), qlm_id( l, m ), qm2m( m ) );
        }
        for( l = 3; l < L; ++l ) {
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] = {};\n", id( l, m ), qlm_id( l, m ), qlm( l, m, ylm( l - 1, m ), ylm( l - 2, m ) ) ) ;
        }
        code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0; c2 = 0; s2 = 0;\n";
        //mについてのループ
        //1. y_m^m (y_m^{-m})を計算
        //2. y_{m+1}^mをy_m^mから計算
        //3. y_l^mをy_{l-1}^m, y_{l-2}^mから計算
        for( m = 1; m < L - 1; ++m ) {
            code << tab() << std::format( "//m = {:03d}\n", m );
            l = m;
            code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
            code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );

            if( m + 1 < L ) {
                l++;
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm1m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }
            if( m + 2 < L ) {
                l++;
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qm2m( m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }
            for( l = m + 3; l < L; ++l ) {
                code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), qlm( l, m, std::format( "qlm[ {:3d} ]", qlm_id( l - 1, m ) ), std::format( "qlm[ {:3d} ]", qlm_id( l - 2, m ) ) ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
                code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );
            }
            //update cosine and sine
            code << tab() << "s2 = s1; c2 = c1; s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
        }
        //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        code << tab() << std::format( "//m = {:03d}\n", m );
        code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
        code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );

        //calculate gradient
        // m = 0 Zonal Harmonics
        code << tab() << "//calculate gradient\n";
        code << tab() << "const glm::vec3 zero = glm::vec3( 0.f, 0.f, 0.f );\n";
        l = m = 0;
        code << tab() << std::format( "glm[ {:3d} ] = zero;\n", id( l, m ) );
        l = 1;
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * ylm[ {0:3d} ];\n", id( l, m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * ylm[ {0:3d} ];\n", id( l, m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * ylm[ {0:3d} ] + {1};\n", id( l, m ), constant( klm( l, m ) ) );
        for( l = 2; l < L; ++l ) {
            code << tab() << std::format( "tmp0 = {:3d} * ylm[ {:3d} ]; tmp1 = {} * qlm[ {:3d} ];\n", l, id( l, m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ) );
            code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp0 - tmp1 );\n", id( l, m ) );
            code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp0 - tmp1 );\n", id( l, m ) );
            code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp0 + {} * qlm[ {:3d} ];\n", id( l, m ), constant( klm( l, m ) / klm( l - 1, m ) * l ), qlm_id( l - 1, m ) );
        }
        code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0;\n";
        for( m = 1; m < L - 1; ++m ) {
            code << tab() << std::format( "//m = {:03d}\n", m );
            l = m;
            if( l == 1 ) {
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ylm[ {:3d} ] - {};\n", id( l, + m ), id( l, + m ), sqrt( 3 / ( 4.0 * std::numbers::pi ) ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ylm[ {:3d} ];\n", id( l, + m ), id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * ylm[ {:3d} ];\n", id( l, + m ), id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ylm[ {:3d} ];\n", id( l, - m ), id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ylm[ {:3d} ] - {};\n", id( l, - m ), id( l, - m ), sqrt( 3 / ( 4.0 * std::numbers::pi ) ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * ylm[ {:3d} ];\n", id( l, - m ), id( l, - m ) );
            } else {
                const double coef = klm( l, l ) / klm( l - 1, l - 1 ) * l * ( 2.0 * l - 1.0 );
                code << tab() << std::format( "tmp0 = {} * ylm[ {:3d} ]; tmp1 = {} * ylm[ {:3d} ]; tmp2 = {} * ylm[ {:3d} ]; tmp3 = {} * ylm[ {:3d} ];\n", l, id( l, + m ), l, id( l, - m ), constant( coef ), id( l - 1, l - 1 ), constant( coef ), id( l - 1, 1 - l ) );
                code << tab() << std::format( "glm[ {:3d} ].x = - x * tmp0 - tmp2;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * tmp0 + tmp3;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp0;\n", id( l, + m ) ); //, constant( klm( l, m ) / klm( l - 1, l ) * ( l + m ) ), id( l - 1, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].x = - x * tmp1 - tmp3;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * tmp1 - tmp2;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp1;\n", id( l, - m ) ); //, constant( klm( l, m ) / klm( l - 1, l ) * ( l + m ) ), id( l - 1, - m ) );
            }

            for( l = m + 1; l < L; ++l ) {
                code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c0; tmp1 = tmp * s0;\n", constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ) );
                code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp2 = tmp * c1; tmp3 = tmp * s1;\n", m, qlm_id( l, m ) );
                code << tab() << std::format( "tmp = {} * ylm[ {:3d} ];\n", l, id( l, + m ) );
                //grad x & y
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp - tmp0 ) + tmp2;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp - tmp0 ) - tmp3;\n", id( l, + m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp + {} * ylm[ {:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, + m ) );

                code << tab() << std::format( "tmp = {} * ylm[ {:3d} ];\n", l, id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].x = - x * ( tmp - tmp1 ) + tmp3;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].y = - y * ( tmp - tmp1 ) + tmp2;\n", id( l, - m ) );
                code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp + {} * ylm[ {:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, - m ) );

                //code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; glm[ {:3d} ].z += tmp * c0; glm[ {:3d} ].z += tmp * s0;\n", constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), qlm_id( l - 1, m ), id( l, + m ), id( l, - m ) );
            }
            code << tab() << "s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
        }
        //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        const double coef = klm( l, l ) / klm( l - 1, l - 1 ) * l * ( 2.0 * l - 1.0 );
        code << tab() << std::format( "tmp0 = {} * ylm[ {:3d} ]; tmp1 = {} * ylm[ {:3d} ]; tmp2 = {} * ylm[ {:3d} ]; tmp3 = {} * ylm[ {:3d} ];\n", l, id( l, + m ), l, id( l, - m ), constant( coef ), id( l - 1, l - 1 ), constant( coef ), id( l - 1, 1 - l ) );
        code << tab() << std::format( "glm[ {:3d} ].x = - x * tmp0 - tmp2;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {:3d} ].y = - y * tmp0 + tmp3;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp0;\n", id( l, + m ) ); //, constant( klm( l, m ) / klm( l - 1, l ) * ( l + m ) ), id( l - 1, + m ) );
        code << tab() << std::format( "glm[ {:3d} ].x = - x * tmp1 - tmp3;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {:3d} ].y = - y * tmp1 - tmp2;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp1;\n", id( l, - m ) ); //, constant( klm( l, m ) / klm( l - 1, l ) * ( l + m ) ), id( l - 1, - m ) );
        /*
        code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c1; tmp1 = tmp * s1; tmp2 = {} * ylm[ {:3d} ]; tmp3 = {} * ylm[ {:3d} ];\n", m, qlm_id( l, m ), l, id( l, + m ), l, id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * tmp2 + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * tmp2 - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * tmp2;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - x * tmp3 + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - y * tmp3 + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - z * tmp3;\n", id( l, - m ) );
        */

        ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //calculate hessian matrix
        //zonal harmonics
        m = 0;
        l = 0;
        code << tab() << std::format( "hlm[ {0:3d} ] = glm::mat3( 0.f );\n", id( l, m ) );
        l = 1;
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - ( ( 1 - x2 ) * ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = - ( - xy * ylm[ {0:3d} ] + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = - ( - xz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].x + x * glm[ {0:3d} ].z );\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 1 ];\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - ( ( 1 - y2 ) * ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = - ( - yz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].y + y * glm[ {0:3d} ].z );\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 2 ];\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 2 ];\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - ( ( 1 - z2 ) * ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, m ) );

        for( l = 2; l < L; ++l ) {
            code << tab() << std::format( "tmp0 = {} * qlm[ {:3d} ]; tmp1 = {} * qlm[ {:3d} ]; tmp2 = {} * qlm[ {:3d} ];\n", constant( klm( l, m ) / klm( l - 1, 1 ) ), qlm_id( l - 1, 1 ), constant( klm( l, m ) / klm( l - 2, 2 ) ), qlm_id( l - 2, 2 ), constant( klm( l, m ) / klm( l - 2, 1 ) * l ), qlm_id( l - 2, 1 ) );
            code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = tmp0 + x2 * tmp1 - {4:3d} * ( ( 1 + {5:3d} * x2 ) * ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, m ), constant( klm( l, m ) ), qlm_id( l - 2, 2 ), qlm_id( l - 1, 1 ), l, l - 2 );
            code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = xy * tmp1 - {3:3d} * ( {4:3d} * xy * ylm[ {0:3d} ] + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, m ), constant( klm( l, m ) ), qlm_id( l - 2, 2 ), l, l - 2 );
            code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = x  * tmp2 - {3:3d} * ( {4:3d} * xz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].x + x * glm[ {0:3d} ].z );\n", id( l, m ), constant( klm( l, m ) * l ), qlm_id( l - 2, 1 ), l, l - 2 );
            code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = tmp0 + y2 * tmp1 - {4:3d} * ( ( 1 + {5:3d} * y2 ) * ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, m ), constant( klm( l, m ) ), qlm_id( l - 2, 2 ), qlm_id( l - 1, 1 ), l, l - 2 );
            code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = y  * tmp2 - {3:3d} * ( {4:3d} * yz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].y + y * glm[ {0:3d} ].z );\n", id( l, m ), constant( klm( l, m ) * l ), qlm_id( l - 2, 1 ), l, l - 2 );
            code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = {1} * qlm[ {2:3d} ] - {3:3d} * ( ( 1 + {4:3d} * z2 ) * ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, m ), constant( klm( l, m ) / klm( l - 2, m ) * l * ( l - 1 ) ), qlm_id( l - 2, 0 ), l, l - 2 );
        }
        //
        code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0; c2 = 0; s2 = 0;\n";
        for( m = 1; m < L - 1; ++m ) {
            //
            for( l = m; l < L; ++l ) {
                ////////////////////////////////////////////////////////////////////////////////////////////////////////
                //set initial values for xx, xy, xz, yy, yz, zz
                code << tab() << std::format( "tmp0 = {0} * qlm[ {1:3d} ]; tmp2 = tmp0 * c2; tmp3 = tmp0 * s2;\n", m * ( m - 1 ), qlm_id( l, m ) );
                code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {0} * ylm[ {2:3d} ];\n", l - 2, id( l, + m ), id( l, - m ) ); //tmp0 = (l-2)ylm, tmp1 = (l-2)yl-m
                //xx
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] =   tmp2 - {1} * ( x2 * tmp0 + ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] =   tmp3 - {1} * ( x2 * tmp1 + ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, - m ), l );
                //xy
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = - tmp3 - {1} * ( xy * tmp0 + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] =   tmp2 - {1} * ( xy * tmp1 + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, - m ), l );
                //yy
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp2 - {1} * ( y2 * tmp0 + ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp3 - {1} * ( y2 * tmp1 + ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, - m ), l );
                //xz
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = - {1} * ( xz * tmp0 + x * glm[ {0:3d} ].z + z * glm[ {0:3d} ].x );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = - {1} * ( xz * tmp1 + x * glm[ {0:3d} ].z + z * glm[ {0:3d} ].x );\n", id( l, - m ), l );
                //yz
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = - {1} * ( yz * tmp0 + y * glm[ {0:3d} ].z + z * glm[ {0:3d} ].y );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = - {1} * ( yz * tmp1 + y * glm[ {0:3d} ].z + z * glm[ {0:3d} ].y );\n", id( l, - m ), l );
                //zz
                code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - {1} * ( z2 * tmp0 + ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, + m ), l );
                code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - {1} * ( z2 * tmp1 + ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, - m ), l );
                // xx, xy, yy
                //Q_{l-1}^{m+1} term
                if( m + 1 <= l - 1 ) {
                    code << tab() << std::format( "tmp  = {} * qlm[ {:3d} ];\n", constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ) );
                    code << tab() << std::format( "tmp0 = tmp * {0} * c1; tmp1 = tmp * {0} * s1; tmp2 = tmp * c0; tmp3 = tmp * s0;\n", m ); //tmp1 = mk_l^mQ_{l-1}^{m+1}c_{m-1}, tmp2 = mk_l^mQ_{l-1}^{m+1}s_{m-1}
                    //xx
                    code << tab() << std::format( "hlm[ {:3d} ][ 0 ][ 0 ] += ( tmp0 + tmp0 ) * x + tmp2;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {:3d} ][ 0 ][ 0 ] += ( tmp1 + tmp1 ) * x + tmp3;\n", id( l, - m ) );
                    //xy
                    code << tab() << std::format( "hlm[ {:3d} ][ 0 ][ 1 ] += tmp0 * y - tmp1 * x;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {:3d} ][ 0 ][ 1 ] += tmp0 * x + tmp1 * y;\n", id( l, - m ) );
                    //yy
                    code << tab() << std::format( "hlm[ {:3d} ][ 1 ][ 1 ] += - ( tmp1 + tmp1 ) * y + tmp2;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {:3d} ][ 1 ][ 1 ] +=   ( tmp0 + tmp0 ) * y + tmp3;\n", id( l, - m ) );
                }
                //Q_{l-2}^{m+2} term
                if( m + 2 <= l - 2 ) {
                    code << tab() << std::format( "tmp = {0} * qlm[ {1:3d} ]; tmp0 = tmp * c0; tmp1 = tmp * s0;\n", constant( klm( l, m ) / klm( l - 2, m + 2 ) ), qlm_id( l - 2, m + 2 ) );
                    //xx
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] += x2 * tmp0;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] += x2 * tmp1;\n", id( l, - m ) );
                    //xy
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] += xy * tmp0;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] += xy * tmp1;\n", id( l, - m ) );
                    //yy
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] += y2 * tmp0;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] += y2 * tmp1;\n", id( l, - m ) );
                }
                //xz, yz
                if( m <= l - 1 && ( ( l + m ) * m ) != 0 ) {
                    code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c1; tmp1 = tmp * s1;\n", constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) * m ), qlm_id( l - 1, m ) );
                    //xz
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] += tmp0;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] += tmp1;\n", id( l, - m ) );
                    //yz
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] -= tmp1;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] += tmp0;\n", id( l, - m ) );
                }
                if( m + 1 <= l - 2 && ( l + m ) != 0 ) {
                    code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; tmp0 = tmp * c0; tmp1 = tmp * s0;\n", constant( klm( l, m ) / klm( l - 2, m + 1 ) * ( l + m ) ), qlm_id( l - 2, m + 1 ) );
                    //xz
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] += tmp0 * x;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] += tmp1 * x;\n", id( l, - m ) );
                    //yz
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] += tmp0 * y;\n", id( l, + m ) );
                    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] += tmp1 * y;\n", id( l, - m ) );
                }
                //zz
                if( m <= l - 2 ) {
                    //if( ( l + m ) * ( l + m - 1 ) != 0 ) code << tab() << std::format( "tmp = {} * qlm[ {:3d} ]; hlm[ {:3d} ][ 2 ][ 2 ] += tmp * c0; hlm[ {:3d} ][ 2 ][ 2 ] += tmp * s0;\n", constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), qlm_id( l - 2, m ), id( l, + m ), id( l, - m ) );
                    if( ( l + m ) * ( l + m - 1 ) != 0 ) code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] += {1} * ylm[ {2:3d} ]; hlm[ {3:3d} ][ 2 ][ 2 ] += {1} * ylm[ {4:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, m ), id( l, - m ), id( l - 2, - m ) );
                }
                //copy symmetric elements yx = xy, zx = xz, zy = yz
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 1 ]; hlm[ {0:3d} ][ 2 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 2 ]; hlm[ {0:3d} ][ 2 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 2 ];\n", id( l, + m ) );
                code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 1 ]; hlm[ {0:3d} ][ 2 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 2 ]; hlm[ {0:3d} ][ 2 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 2 ];\n", id( l, - m ) );
            }
            //update cosine and sine
            code << tab() << std::format( "s2 = s1; c2 = c1; s1 = s0; c1 = c0;\n" );
            code << tab() << "c0 = " << update_cos( "c1", "s1" ) << ";\n";
            code << tab() << "s0 = " << update_sin( "c1", "s1" ) << ";\n";
        }
        //last pair hlm( L - 1, L - 1 ), hlm( L - 1, - ( L - 1 ) )
        l = L - 1;
        m = l;
        code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {0} * ylm[ {2:3d} ];\n", constant( klm( l, m ) / klm( l - 2, m - 2 ) * l * ( l - 1 ) * ( 1.0 - 2.0 * l ) * ( 3.0 - 2.0 * l ) ), id( l - 2, m - 2 ), id( l - 2, 2 - m ) );
        //code << tab() << std::format( "tmp0 = {0} * qlm[ {1:3d} ] * c2; tmp1 = {0} * qlm[ {1:3d} ] * s2;\n", m * ( m - 1 ), qlm_id( l, m ) );
        //xx
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = tmp0 - {3:3d} * ( ( 1 + {4:3d} * x2 ) * ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, + m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = tmp1 - {3:3d} * ( ( 1 + {4:3d} * x2 ) * ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, - m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        //xy
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp1 - {3:3d} * ( {4:3d} * xy * ylm[ {0:3d} ] + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, + m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] =   tmp0 - {3:3d} * ( {4:3d} * xy * ylm[ {0:3d} ] + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, - m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        //xz
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - {1:3d} * ( {2:3d} * xz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].x + x * glm[ {0:3d} ].z );\n", id( l, + m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - {1:3d} * ( {2:3d} * xz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].x + x * glm[ {0:3d} ].z );\n", id( l, - m ), l, l - 2 );
        //yy
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp0 - {3:3d} * ( ( 1 + {4:3d} * y2 ) * ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, + m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp1 - {3:3d} * ( ( 1 + {4:3d} * y2 ) * ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, - m ), m * ( m - 1 ), qlm_id( l, m ), l, l - 2 );
        //yz
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - {1:3d} * ( {2:3d} * yz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].y + y * glm[ {0:3d} ].z );\n", id( l, + m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - {1:3d} * ( {2:3d} * yz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].y + y * glm[ {0:3d} ].z );\n", id( l, - m ), l, l - 2 );
        //zz
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - {1:3d} * ( ( 1 + {2:3d} * z2 ) * ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, + m ), l, l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - {1:3d} * ( ( 1 + {2:3d} * z2 ) * ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, - m ), l, l - 2 );
        code << "}\n";
        return code.str();
    }
#pragma endregion

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Spherical Harmonicsのhessianを計算するコードの生成 乗算数を減らすように工夫
    std::string generate_hessian( const std::string &function_name = "shhessian" ) const;

};

#include"stable/shh_code_generator-impl.hpp"