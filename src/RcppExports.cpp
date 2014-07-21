// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// gv_abcd2vis
ComplexVector gv_abcd2vis(NumericVector ii, NumericVector kx = NumericVector::create(4), NumericVector sx = NumericVector::create(4), NumericVector fx = NumericVector::create(1));
RcppExport SEXP gRavity_gv_abcd2vis(SEXP iiSEXP, SEXP kxSEXP, SEXP sxSEXP, SEXP fxSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type ii(iiSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type kx(kxSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type sx(sxSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type fx(fxSEXP );
        ComplexVector __result = gv_abcd2vis(ii, kx, sx, fx);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_const
List gv_const();
RcppExport SEXP gRavity_gv_const() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = gv_const();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_readfits
NumericVector gv_readfits(CharacterVector fits_name, CharacterVector hdu_name = "", CharacterVector col_name = "");
RcppExport SEXP gRavity_gv_readfits(SEXP fits_nameSEXP, SEXP hdu_nameSEXP, SEXP col_nameSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< CharacterVector >::type fits_name(fits_nameSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type hdu_name(hdu_nameSEXP );
        Rcpp::traits::input_parameter< CharacterVector >::type col_name(col_nameSEXP );
        NumericVector __result = gv_readfits(fits_name, hdu_name, col_name);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_sci_img2pix
NumericVector gv_sci_img2pix(NumericMatrix img, NumericVector idx, IntegerVector cnr, IntegerVector isYjunc = 0, IntegerVector n_wd = 3);
RcppExport SEXP gRavity_gv_sci_img2pix(SEXP imgSEXP, SEXP idxSEXP, SEXP cnrSEXP, SEXP isYjuncSEXP, SEXP n_wdSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type img(imgSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type idx(idxSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type cnr(cnrSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type isYjunc(isYjuncSEXP );
        Rcpp::traits::input_parameter< IntegerVector >::type n_wd(n_wdSEXP );
        NumericVector __result = gv_sci_img2pix(img, idx, cnr, isYjunc, n_wd);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_sci_pix2vis_unsorted
List gv_sci_pix2vis_unsorted(NumericVector pixels, NumericVector rdnoiz, List v2pms);
RcppExport SEXP gRavity_gv_sci_pix2vis_unsorted(SEXP pixelsSEXP, SEXP rdnoizSEXP, SEXP v2pmsSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericVector >::type pixels(pixelsSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rdnoiz(rdnoizSEXP );
        Rcpp::traits::input_parameter< List >::type v2pms(v2pmsSEXP );
        List __result = gv_sci_pix2vis_unsorted(pixels, rdnoiz, v2pms);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_solvels
List gv_solvels(NumericMatrix rA, NumericVector rb, NumericVector rw);
RcppExport SEXP gRavity_gv_solvels(SEXP rASEXP, SEXP rbSEXP, SEXP rwSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< NumericMatrix >::type rA(rASEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rb(rbSEXP );
        Rcpp::traits::input_parameter< NumericVector >::type rw(rwSEXP );
        List __result = gv_solvels(rA, rb, rw);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// gv_vis2gd
List gv_vis2gd(ComplexVector vis);
RcppExport SEXP gRavity_gv_vis2gd(SEXP visSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< ComplexVector >::type vis(visSEXP );
        List __result = gv_vis2gd(vis);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP gRavity_rcpp_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = rcpp_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// timesTwo
NumericVector timesTwo();
RcppExport SEXP gRavity_timesTwo() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        NumericVector __result = timesTwo();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
