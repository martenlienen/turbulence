#include <gtest/gtest.h>
#include "stencils/StencilFunctions.h"

TEST(X, x) { ASSERT_EQ(5, 2 + 3); }

TEST(stencil, FdnududxdxTest) {
  FLOAT lv[27 * 3];
  FLOAT lm[27 * 3];
  FLOAT ln[27 * 3];

  lv[mapd(-1, +0, +0, +0)] = 1.0;
  lv[mapd(+0, +0, +0, +0)] = 2.0;
  lv[mapd(+1, +0, +0, +0)] = 3.0;

  lm[mapd(+0, +0, +0, +0)] = 1.0;
  lm[mapd(+1, +0, +0, +0)] = 1.0;

  ln[mapd(+0, +0, +0, +0)] = 1.0;
  ln[mapd(+1, +0, +0, +0)] = 1.0;

  ASSERT_EQ(+0.0, Fdnududxdx(lv, lm, ln));

  lv[mapd(+1, +0, +0, +0)] = 1.0;
  ASSERT_EQ(-2.0, Fdnududxdx(lv, lm, ln));
}

TEST(stencil, FdnududydyTest) {
  FLOAT lv[27 * 3];
  FLOAT lm[27 * 3];
  FLOAT ln[27 * 3];

  lv[mapd(+0, -1, +0, +U)] = 1.0;
  lv[mapd(+0, +0, +0, +U)] = 2.0;
  lv[mapd(+0, +1, +0, +U)] = 3.0;

  lm[mapd(+0, -1, +0, +X)] = 1.0;
  lm[mapd(+1, -1, +0, +X)] = 1.0;
  lm[mapd(+0, +0, +0, +X)] = 1.0;
  lm[mapd(+1, +0, +0, +X)] = 1.0;
  lm[mapd(+0, +1, +0, +X)] = 1.0;
  lm[mapd(+1, +1, +0, +X)] = 1.0;
  lm[mapd(+0, -1, +0, +Y)] = 1.0;
  lm[mapd(+1, -1, +0, +Y)] = 1.0;
  lm[mapd(+0, +0, +0, +Y)] = 1.0;
  lm[mapd(+1, +0, +0, +Y)] = 1.0;
  lm[mapd(+0, +1, +0, +Y)] = 1.0;
  lm[mapd(+1, +1, +0, +Y)] = 1.0;

  ln[mapd(+0, -1, +0, +X)] = 1.0;
  ln[mapd(+1, -1, +0, +X)] = 1.0;
  ln[mapd(+0, +0, +0, +X)] = 1.0;
  ln[mapd(+1, +0, +0, +X)] = 1.0;
  ln[mapd(+0, +1, +0, +X)] = 1.0;
  ln[mapd(+1, +1, +0, +X)] = 1.0;

  ASSERT_EQ(+0.0, Fdnududydy(lv, lm, ln));

  lv[mapd(+0, +1, +0, +0)] = 1.0;
  ASSERT_EQ(-2.0, Fdnududydy(lv, lm, ln));
}

TEST(stencil, FdnududydxTest) {
  FLOAT lv[27 * 3];
  FLOAT lm[27 * 3];
  FLOAT ln[27 * 3];

  lv[mapd(+0, -1, +0, +U)] = 1.0;
  lv[mapd(+0, +0, +0, +U)] = 2.0;
  lv[mapd(+0, +1, +0, +U)] = 3.0;

  lm[mapd(+0, -1, +0, +X)] = 1.0;
  lm[mapd(+1, -1, +0, +X)] = 1.0;
  lm[mapd(+0, +0, +0, +X)] = 1.0;
  lm[mapd(+1, +0, +0, +X)] = 1.0;
  lm[mapd(+0, +1, +0, +X)] = 1.0;
  lm[mapd(+1, +1, +0, +X)] = 1.0;
  lm[mapd(+0, -1, +0, +Y)] = 1.0;
  lm[mapd(+1, -1, +0, +Y)] = 1.0;
  lm[mapd(+0, +0, +0, +Y)] = 1.0;
  lm[mapd(+1, +0, +0, +Y)] = 1.0;
  lm[mapd(+0, +1, +0, +Y)] = 1.0;
  lm[mapd(+1, +1, +0, +Y)] = 1.0;

  ln[mapd(+0, -1, +0, +X)] = 1.0;
  ln[mapd(+1, -1, +0, +X)] = 1.0;
  ln[mapd(+0, +0, +0, +X)] = 1.0;
  ln[mapd(+1, +0, +0, +X)] = 1.0;
  ln[mapd(+0, +1, +0, +X)] = 1.0;
  ln[mapd(+1, +1, +0, +X)] = 1.0;

  ASSERT_EQ(+0.0, Fdnududydy(lv, lm, ln));

  lv[mapd(+0, +1, +0, +0)] = 1.0;
  ASSERT_EQ(-2.0, Fdnududydy(lv, lm, ln));
}

TEST(stencil, biLinearInterpolationTest) {
  ASSERT_EQ(1.0, biLinearInterpolation(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0));

  // CORNERS
  ASSERT_EQ(1.0, biLinearInterpolation(0.0, 1.0, 0.0, 1.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(2.0, biLinearInterpolation(1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(3.0, biLinearInterpolation(1.0, 0.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(4.0, biLinearInterpolation(0.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0));

  // CENTER OF EDGES
  ASSERT_EQ(1.5, biLinearInterpolation(1.0, 1.0, 0.0, 1.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(2.5, biLinearInterpolation(1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(3.5, biLinearInterpolation(1.0, 1.0, 1.0, 0.0, 1.0, 2.0, 3.0, 4.0));

  ASSERT_EQ(2.5, biLinearInterpolation(0.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0));

  // CENTRE
  ASSERT_EQ(2.5, biLinearInterpolation(1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 3.0, 4.0));
}