#pragma once
// util_stable.hpp
//コード生成器で使用するユーティリティ関数
#include<iostream>
#include<format>
#include<print>
#include<numbers>

namespace util {

//定数を文字列に変換(小数点以下16桁)
inline std::string constant( const double x )
{
    return std::format( "{:.16f}f", x );
}


//代入文
inline std::string assign( const std::string &lhs, const std::string &rhs )
{
    return std::format( "{} = {}", lhs, rhs );
}

//演算 opは演算子(+,-,*,%)
inline std::string op( const std::string &op0, const std::string &op1, const std::string &op )
{
    return std::format( "{} {} {}", op0, op, op1 );
}

//加算
inline std::string add( const std::string &op0, const std::string &op1 )
{
    return op( op0, op1, "+" );
}

//減算
inline std::string sub( const std::string &op0, const std::string &op1 )
{
    return op( op0, op1, "-" );
}

//乗算
inline std::string mul( const std::string &op0, const std::string &op1 )
{
    return op( op0, op1, "*" );
}

//配列変数"variable[ id ]"の文字列を返す
inline std::string array( const std::string &variable, const int id )
{
    return std::format( "{}[ {:3} ]", variable, id );
}

//体球調和関数の配列"ylm[ id ]"の文字列を返す
inline std::string ylm( const int l, const int m )
{
    return array( "ylm", l * ( l + 1 ) + m );
}

//体球調和関数の勾配の配列"glm[ id ]"の文字列を返す
inline std::string glm( const int l, const int m )
{
    return array( "glm", l * ( l + 1 ) + m );
}

//タブ
inline std::string tab()
{
    return "    ";
}

//カッコ()を文字列strの前後につけた文字列( str )を返す
inline std::string parenthesis( const std::string &str )
{
    return std::format( "( {} )", str );
}

//SHの正規化定数K_l^mを返す
inline double Klm( const int l, const int m )
{
    const unsigned int absm = ( m >= 0 ) ? m : - m;
    double v = 1.0;
    for( unsigned int k = l + absm; k > ( l - absm ); k-- ) v *= k;
    return sqrt( ( 2.0 * l + 1.0 ) / ( 4 * std::numbers::pi * v ) );
}

inline double Qmm( const int m )
{
    double val = 1.0;
    for( int n {}; n <= m; ++n ) val *= ( 1.0 - 2.0 * n );
    return val;
}

//SHの正規化定数klm = \sqrt(2)Klm (m>0) : Klm (m=0)を返す
inline double klm( const int l, const int m )
{
    return ( m == 0 ) ? Klm( l, 0 ) : sqrt( 2.0 ) * Klm( l, m );
}

//qmm is constant klm * Qmm
inline double qmm( const int m )
{
    return klm( m, m ) * Qmm( m );
}

//q_{m+1}^m = k_{m+1}^m/k_m^m * (2m+1)*z*qmm
inline std::string qm1m( const int m )
{
    return std::format( "{} * z", constant( klm( m + 1, m ) / klm( m, m ) * ( 2 * m + 1 ) * qmm( m ) ) );
}

//inline std::string Qm1m( const int m )
//{
//    return std::format( "{} * z", constant( ( 2 * m + 1 ) * Qmm( m ) ) );
//}

//q_{m+2}^m = k_{m+2}^m/k_{m}^m(2m+1)qmm((2m+3)z^2 - r^2)/2
inline std::string qm2m( const int m )
{
    const std::string coef0 = constant( klm( m + 2, m ) / klm( m, m ) * ( 2 * m + 1 ) * ( 2 * m + 3 ) * qmm( m ) / 2.0 );
    const auto val = klm( m + 2, m ) / klm( m, m ) * ( 2 * m + 1 ) * qmm( m ) / 2.0;
    if( val < 0.0 ) {
        return std::format( "{} * z2 + {}", coef0, constant( - val ) );
    } else {
        return std::format( "{} * z2 - {}", coef0, constant( val ) );
    }
}

inline std::string Qm2m( const int m )
{
    const std::string coef0 = constant( ( 2 * m + 1 ) * ( 2 * m + 3 ) * Qmm( m ) / 2.0 );
    const auto val = ( 2 * m + 1 ) * Qmm( m ) / 2.0;
    if( val < 0.0 ) {
        return std::format( "{} * z2 + {}", coef0, constant( - val ) );
    } else {
        return std::format( "{} * z2 - {}", coef0, constant( val ) );
    }
}

//q_l^m = (2l-1)*z*q_{l-1}^m * k_l^m/k_{l-1}^m / ( l - m ) - ( l + m - 1 ) * q_{l-2}^m * k_l^m/k_{l-2}^m / ( l - m )
inline std::string qlm( const int l, const int m, const std::string &ql1m, const std::string &ql2m )
{
    const std::string coef0 = constant( ( 2 * l - 1.0 ) / ( l - m ) * klm( l, m ) / klm( l - 1, m ) );
    const std::string coef1 = constant( ( l + m - 1.0 ) / ( l - m ) * klm( l, m ) / klm( l - 2, m ) );
    return std::format( "{} * z * {} - {} * {}", coef0, ql1m, coef1, ql2m );
}

inline std::string Qlm( const int l, const int m, const std::string &ql1m, const std::string &ql2m )
{
    const std::string coef0 = constant( ( 2 * l - 1.0 ) / ( l - m ) );
    const std::string coef1 = constant( ( l + m - 1.0 ) / ( l - m ) );
    return std::format( "{} * z * {} - {} * {}", coef0, ql1m, coef1, ql2m );
}

// cm = x * c{m-1} - y * s{m-1}
inline std::string update_cos( const std::string &sc, const std::string &ss )
{
    return std::format( "x * {} - y * {}", sc, ss );
}

// sm = y * c{m-1} + x * s{m-1}
inline std::string update_sin( const std::string &sc, const std::string &ss )
{
    return std::format( "y * {} + x * {}", sc, ss );
}


}

