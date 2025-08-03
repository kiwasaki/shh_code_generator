
template< size_t L > std::string shh_code_generator< L >::generate_hessian( const std::string &function_name ) const
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
    code << tab() << "float c0, c1, cm, cs, s0, s1, sm, ss, tmp, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, lx, ly, lz;\n";
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
    code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0;\n";
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
        code << tab() << "s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n";
    }
    //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
    l = L - 1;
    m = l;
    code << tab() << std::format( "//m = {:03d}\n", m );
    code << tab() << std::format( "qlm[ {:3d} ] = {};\n", qlm_id( l, m ), constant( qmm( m ) ) );
    code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * c0;\n", id( l, + m ), qlm_id( l, m ) );
    code << tab() << std::format( "ylm[ {:3d} ] = qlm[ {:3d} ] * s0;\n", id( l, - m ), qlm_id( l, m ) );

    //calculate gradient & hessian
    // m = 0 Zonal Harmonics
    code << tab() << "//calculate gradient\n";
    //code << tab() << "const glm::vec3 zero = glm::vec3( 0.f, 0.f, 0.f );\n";
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //(0,0)
    l = m = 0;
    code << tab() << std::format( "glm[ {:3d} ] = {{}};\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {:3d} ] = glm::mat3( 0.f );\n", id( l, m ) );

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //(1,0)
    l = 1;
    code << tab() << std::format( "lx = x; ly = y; lz = z;\n" );
    code << tab() << std::format( "glm[ {0:3d} ].x = - x * ylm[ {0:3d} ];\n", id( l, m ) );
    code << tab() << std::format( "glm[ {0:3d} ].y = - y * ylm[ {0:3d} ];\n", id( l, m ) );
    code << tab() << std::format( "glm[ {0:3d} ].z = - z * ylm[ {0:3d} ] + {1};\n", id( l, m ), constant( klm( l, m ) ) );

    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - ( ( 1 - x2 ) * ylm[ {0:3d} ] + 2 * x * glm[ {0:3d} ].x );\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = - ( - xy * ylm[ {0:3d} ] + y * glm[ {0:3d} ].x + x * glm[ {0:3d} ].y );\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = - ( - xz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].x + x * glm[ {0:3d} ].z );\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 1 ];\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - ( ( 1 - y2 ) * ylm[ {0:3d} ] + 2 * y * glm[ {0:3d} ].y );\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = - ( - yz * ylm[ {0:3d} ] + z * glm[ {0:3d} ].y + y * glm[ {0:3d} ].z );\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 0 ] = hlm[ {0:3d} ][ 0 ][ 2 ];\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 2 ];\n", id( l, m ) );
    code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - ( ( 1 - z2 ) * ylm[ {0:3d} ] + 2 * z * glm[ {0:3d} ].z );\n", id( l, m ) );
    code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // (2,0), (3,0) ...(L-1,0)
    for( l = 2; l < L; ++l ) {
        //gradient
        code << tab() << std::format( "tmp = {:3d} * ylm[ {:3d} ]; tmp0 = {} * qlm[ {:3d} ] - tmp;\n", l, id( l, + m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m  + 1 ) );
        code << tab() << std::format( "glm[ {:3d} ].x = x * tmp0;\n", id( l, m ) );
        code << tab() << std::format( "glm[ {:3d} ].y = y * tmp0;\n", id( l, m ) );
        code << tab() << std::format( "glm[ {:3d} ].z = - z * tmp + {} * qlm[ {:3d} ];\n", id( l, m ), constant( klm( l, m ) / klm( l - 1, m ) * l ), qlm_id( l - 1, m ) );
        //hessian
        code << tab() << std::format( "tmp1 = {} * tmp; tmp2 = {} * qlm[ {:3d} ] - tmp1; tmp3 = {} * qlm[ {:3d} ] - glm[ {:3d} ].z;\n", l - 2, constant( klm( l, m ) / klm( l - 2, m + 2 ) ), qlm_id( l - 2, m + 2 ), constant( klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ), id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = x2 * tmp2 + tmp0 - 2 * lx * glm[ {0:3d} ].x;\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = xy * tmp2 - ly * glm[ {0:3d} ].x - lx * glm[ {0:3d} ].y;\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = y2 * tmp2 + tmp0 - 2 * ly * glm[ {0:3d} ].y;\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - xz * tmp1 + lx * tmp3 - lz * glm[ {0:3d} ].x;\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - yz * tmp1 + ly * tmp3 - lz * glm[ {0:3d} ].y;\n", id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - z2 * tmp1 - 2 * lz * glm[ {0:3d} ].z + {1} * qlm[ {2:3d} ] - tmp;\n", id( l, m ), constant( klm( l, m ) / klm( l - 2, m ) * l * ( l - 1 ) ), qlm_id( l - 2, m ) );
        code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //(l,m)
    // l = m
    const auto block_m0 = [&]() {
        //gradient
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cm; tmp1 = qlm[ {0:3d} ] * sm;\n", qlm_id( l, m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - lx * ylm[ {0:3d} ] + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - ly * ylm[ {0:3d} ] - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ]       ;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - lx * ylm[ {0:3d} ] + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - ly * ylm[ {0:3d} ] + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ]       ;\n", id( l, - m ) );
        //hessian tmp0 = l*ylm, tmp1 = l(l-2)ylm, tmp2 = qlm*cm'', tmp3 = qlm*sm''
        code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {2} * tmp0; tmp2 = qlm[ {3:3d} ] * cs; tmp3 = qlm[ {3:3d} ] * ss;\n", l, id( l, + m ), l - 2, qlm_id( l, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp1 * x2 - tmp0 + tmp2 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp1 * xy        - tmp3 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp1 * y2 - tmp0 - tmp2 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp1 * xz               - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp1 * yz               - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp1 * z2 - tmp0        - ( lz + lz ) * glm[ {0:3d} ].z;\n", id( l, + m ) );
        code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {2} * tmp0;\n", l, id( l, - m ), l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp1 * x2 - tmp0 + tmp3 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp1 * xy        + tmp2 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp1 * y2 - tmp0 - tmp3 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp1 * xz               - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp1 * yz               - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp1 * z2 - tmp0        - ( lz + lz ) * glm[ {0:3d} ].z;\n", id( l, - m ) );
        //update lx, ly, lz
        if( l < L - 1 ) code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    };

    // l = m + 1
    const auto block_m1 = [&]() {
        //gradient
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cm; tmp1 = qlm[ {0:3d} ] * sm; tmp2 = {1} * ylm[ {2:3d} ]; tmp3 = {1} * ylm[ {3:3d} ];\n", qlm_id( l, m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, + m ), id( l - 1, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - lx * ylm[ {0:3d} ] + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - ly * ylm[ {0:3d} ] - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + tmp2;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = - lx * ylm[ {0:3d} ] + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = - ly * ylm[ {0:3d} ] + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + tmp3;\n", id( l, - m ) );
        //hessian tmp0 = l*ylm, tmp1 = l(l-2)ylm, tmp2 = qlm*cm'', tmp3 = qlm*sm''
        code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {2} * tmp0; tmp2 = qlm[ {3:3d} ] * cs; tmp3 = qlm[ {3:3d} ] * ss; tmp4 = {4} * qlm[ {5:3d} ]; tmp5 = tmp4 * sm; tmp4 *= cm;\n", l, id( l, + m ), l - 2, qlm_id( l, m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), qlm_id( l - 1, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp1 * x2 - tmp0 + tmp2 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp1 * xy        - tmp3 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp1 * y2 - tmp0 - tmp2 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp1 * xz        + tmp4 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp1 * yz        - tmp5 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp1 * z2 - tmp0        - ( lz + lz ) * glm[ {0:3d} ].z;\n", id( l, + m ) );
        code << tab() << std::format( "tmp0 = {0} * ylm[ {1:3d} ]; tmp1 = {2} * tmp0;\n", l, id( l, - m ), l - 2 );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp1 * x2 - tmp0 + tmp3 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp1 * xy        + tmp2 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp1 * y2 - tmp0 - tmp3 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp1 * xz        + tmp5 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp1 * yz        + tmp4 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp1 * z2 - tmp0        - ( lz + lz ) * glm[ {0:3d} ].z;\n", id( l, - m ) );
        //update lx, ly, lz
        if( l < L - 1 ) code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    };
    //l = m + 2
    const auto block_m2 = [&]() {
        //gradient
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cm; tmp1 = qlm[ {0:3d} ] * sm; tmp = {1} * qlm[ {2:3d} ]; tmp2 = tmp * c0 - {3} * ylm[ {4:3d} ]; tmp3 = tmp * s0 - {3} * ylm[ {5:3d} ];\n", qlm_id( l, m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ), l, id( l, + m ), id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp2 * x + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp2 * y - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp3 * x + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp3 * y + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, - m ) );
        //hessian
        // tmp0 = qlm*cm'', tmp1 = qlm*sm'', tmp2 = tmp*cm', tmp3 = tmp*sm'
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cs; tmp1 = qlm[ {0:3d} ] * ss; tmp4 = tmp * cm; tmp5 = tmp * sm; tmp = {1} * ylm[ {2:3d} ]; tmp6 = {3} * qlm[ {4:3d} ]; tmp7 = tmp6 * sm; tmp6 *= cm;\n", qlm_id( l, m ), l * ( l - 2 ), id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), qlm_id( l - 1, m )  );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp * x2 + ( tmp4 + tmp4 ) * x + tmp0 + tmp2 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp * xy - tmp5 * x + tmp4 * y - tmp1 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp * y2 - ( tmp5 + tmp5 ) * y - tmp0 + tmp2 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp * xz        + tmp6 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp * yz        - tmp7 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, + m ), l );
        code << tab() << std::format( "tmp = {0} * ylm[ {1:3d} ];\n", l * ( l - 2 ), id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp * x2 + ( tmp5 + tmp5 ) * x + tmp1 + tmp3 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp * xy + tmp4 * x + tmp5 * y + tmp0 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp * y2 + ( tmp4 + tmp4 ) * y - tmp1 + tmp3 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp * xz        + tmp7 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp * yz        + tmp6 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, - m ), l );
        //update lx, ly, lz
        if( l < L - 1 ) code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    };
    // l = m + 3
    const auto block_m3 = [&]() {
        //gradient
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cm; tmp1 = qlm[ {0:3d} ] * sm; tmp = {1} * qlm[ {2:3d} ]; tmp2 = tmp * c0 - {3} * ylm[ {4:3d} ]; tmp3 = tmp * s0 - {3} * ylm[ {5:3d} ];\n", qlm_id( l, m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ), l, id( l, + m ), id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp2 * x + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp2 * y - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp3 * x + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp3 * y + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, - m ) );
        //hessian
        // tmp0 = qlm*cm'', tmp1 = qlm*sm'', tmp2 = tmp*cm', tmp3 = tmp*sm'
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cs; tmp1 = qlm[ {0:3d} ] * ss; tmp4 = tmp * cm; tmp5 = tmp * sm; tmp = {1} * ylm[ {2:3d} ]; tmp6 = {3} * qlm[ {4:3d} ]; tmp7 = tmp6 * sm; tmp6 *= cm;\n", qlm_id( l, m ), l * ( l - 2 ), id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), qlm_id( l - 1, m )  );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp * x2 + ( tmp4 + tmp4 ) * x + tmp0 + tmp2 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp * xy - tmp5 * x + tmp4 * y - tmp1 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp * y2 - ( tmp5 + tmp5 ) * y - tmp0 + tmp2 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "tmp2 = {} * qlm[ {:3d} ] * c0;\n", constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp * xz + tmp2 * x + tmp6 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, + m ), constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp * yz + tmp2 * y - tmp7 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, + m ), constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, + m ), l );
        code << tab() << std::format( "tmp = {0} * ylm[ {1:3d} ]; tmp2 = {2} * qlm[ {3:3d} ] * s0;\n", l * ( l - 2 ), id( l, - m ), constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = - tmp * x2 + ( tmp5 + tmp5 ) * x + tmp1 + tmp3 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = - tmp * xy + tmp4 * x + tmp5 * y + tmp0 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = - tmp * y2 + ( tmp4 + tmp4 ) * y - tmp1 + tmp3 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp * xz + tmp2 * x + tmp7 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp * yz + tmp2 * y + tmp6 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, - m ), l );
        //update lx, ly, lz
        if( l < L - 1 ) code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    };
    // l >= m + 4
    const auto block = [&]() {
        //gradient
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cm; tmp1 = qlm[ {0:3d} ] * sm; tmp = {1} * qlm[ {2:3d} ]; tmp2 = tmp * c0 - {3} * ylm[ {4:3d} ]; tmp3 = tmp * s0 - {3} * ylm[ {5:3d} ];\n", qlm_id( l, m ), constant( klm( l, m ) / klm( l - 1, m + 1 ) ), qlm_id( l - 1, m + 1 ), l, id( l, + m ), id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp2 * x + tmp0;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp2 * y - tmp1;\n", id( l, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, + m ) );
        code << tab() << std::format( "glm[ {0:3d} ].x = tmp3 * x + tmp1;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].y = tmp3 * y + tmp0;\n", id( l, - m ) );
        code << tab() << std::format( "glm[ {0:3d} ].z = - lz * ylm[ {0:3d} ] + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 1, m ) * ( l + m ) ), id( l - 1, - m ) );
        //hessian
        // tmp0 = qlm*cm'', tmp1 = qlm*sm'', tmp2 = tmp*cm', tmp3 = tmp*sm'
        //code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cs; tmp1 = qlm[ {0:3d} ] * ss; tmp4 = tmp * cm; tmp5 = tmp * sm; tmp = {3} * qlm[ {4:3d} ] * c0 - {1} * ylm[ {2:3d} ]; tmp7 = tmp6 * sm; tmp6 *= cm;\n", qlm_id( l, m ), l * ( l - 2 ), id( l, + m ), constant( klm( l, m ) / klm( l - 2, m + 2 ) ), qlm_id( l - 2, m + 2 ) );
        code << tab() << std::format( "tmp0 = qlm[ {0:3d} ] * cs; tmp1 = qlm[ {0:3d} ] * ss; tmp4 = tmp * cm; tmp5 = tmp * sm; tmp = {3} * qlm[ {4:3d} ] * c0 - {1} * ylm[ {2:3d} ]; tmp6 = {5} * qlm[ {6:3d} ]; tmp7 = tmp6 * sm; tmp6 *= cm;\n", qlm_id( l, m ), l * ( l - 2 ), id( l, m ), constant( klm( l, m ) / klm( l - 2, m + 2 ) ), qlm_id( l - 2, m + 2 ), constant( ( l + m ) * klm( l, m ) / klm( l - 1, m ) ), qlm_id( l - 1, m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = tmp * x2 + ( tmp4 + tmp4 ) * x + tmp0 + tmp2 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = tmp * xy - tmp5 * x + tmp4 * y - tmp1 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = tmp * y2 - ( tmp5 + tmp5 ) * y - tmp0 + tmp2 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, + m ) );

        code << tab() << std::format( "tmp = {} * ylm[ {:3d} ]; tmp2 = {} * qlm[ {:3d} ] * c0;\n", l * ( l - 2 ), id( l, + m ), constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp * xz + tmp2 * x + tmp6 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp * yz + tmp2 * y - tmp7 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, + m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, + m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, + m ), l );

        code << tab() << std::format( "tmp = {2} * qlm[ {3:3d} ] * s0 - {0} * ylm[ {1:3d} ]; tmp2 = {4} * ylm[ {5:3d} ];\n", l * ( l - 2 ), id( l, - m ), constant( klm( l, m ) / klm( l - 2, m + 2 ) ), qlm_id( l - 2, m + 2 ), l * ( l - 2 ), id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 0 ] = tmp * x2 + ( tmp5 + tmp5 ) * x + tmp1 + tmp3 - ( lx + lx ) * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 1 ] = hlm[ {0:3d} ][ 1 ][ 0 ] = tmp * xy + tmp4 * x + tmp5 * y + tmp0 - lx * glm[ {0:3d} ].y - ly * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 1 ] = tmp * y2 + ( tmp4 + tmp4 ) * y - tmp1 + tmp3 - ( ly + ly ) * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "tmp = {} * qlm[ {:3d} ] * s0;\n", constant( ( l + m ) * klm( l, m ) / klm( l - 2, m + 1 ) ), qlm_id( l - 2, m + 1 ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 0 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 0 ] = - tmp2 * xz + tmp * x + tmp7 - lx * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].x;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 1 ][ 2 ] = hlm[ {0:3d} ][ 2 ][ 1 ] = - tmp2 * yz + tmp * y + tmp6 - ly * glm[ {0:3d} ].z - lz * glm[ {0:3d} ].y;\n", id( l, - m ) );
        code << tab() << std::format( "hlm[ {0:3d} ][ 2 ][ 2 ] = - tmp2 * z2 - {3} * ylm[ {0:3d} ] - ( lz + lz ) * glm[ {0:3d} ].z + {1} * ylm[ {2:3d} ];\n", id( l, - m ), constant( klm( l, m ) / klm( l - 2, m ) * ( l + m ) * ( l + m - 1 ) ), id( l - 2, - m ), l );
        //update lx, ly, lz
        if( l < L - 1 ) code << tab() << std::format( "lx += x; ly += y; lz += z;\n" );
    };

    code << tab() << "c0 = x; s0 = y; c1 = 1; s1 = 0; cm = 1; sm = 0; cs = 0; ss = 0;\n";
    for( m = 1; m < L - 1; ++m ) {
        code << tab() << std::format( "//m = {:03d}\n", m );
        code << tab() << std::format( "lx = {0} * x; ly = {0} * y; lz = {0} * z;\n", m );
        // l = m
        l = m    ; if( l < L ) block_m0(); //{ block_m0(); //std::println( "{}", id( l, m ) ); }
        l = m + 1; if( l < L ) block_m1(); //{ block_m1(); std::println( "{}", id( l, m ) ); }
        l = m + 2; if( l < L ) block_m2(); //{ block_m2(); std::println( "{}", id( l, m ) ); }
        l = m + 3; if( l < L ) block_m3(); //{ block_m3(); std::println( "{}", id( l, m ) ); }
        for( l = m + 4; l < L; ++l ) block();
        //update c0, s0, c1, s1, c2, s2, cm, sm
        code << tab() << std::format( "cs = {0} * c1; ss = {0} * s1; cm = {1} * c0; sm = {1} * s0; s1 = s0; c1 = c0; c0 = x * c1 - y * s1; s0 = y * c1 + x * s1;\n", ( m + 1 ) * m, m + 1 );
    }
    //last pair sh( L - 1, L - 1 ), sh( L - 1, - ( L - 1 ) )
    l = m = L - 1;
    code << tab() << std::format( "lx = {0} * x; ly = {0} * y; lz = {0} * z;\n", m );
    block_m0();
#if 0
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
#endif
#if 0
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
#endif
    code << "}\n";
    return code.str();
}

