#include <gtest/gtest.h>
#include "clipper2/clipper.h"
#include "ClipFileLoad.h"

inline Clipper2Lib::PathD MakeRandomPath(int width, int height, unsigned vertCnt)
{
  Clipper2Lib::PathD result;
  result.reserve(vertCnt);
  for (unsigned i = 0; i < vertCnt; ++i)
    result.push_back(Clipper2Lib::PointD(double(rand()) / RAND_MAX * width, double(rand()) / RAND_MAX * height));
  return result;
}

template <size_t N>
inline bool IsInList(int num, const int (&intArray)[N])
{  
  const int* list = &intArray[0];
  for (int cnt = N; cnt; --cnt)
    if (num == *list++) return true;
  return false;
}

TEST(Clipper2Tests, TestMultiplePolygons)
{
  std::ifstream ifs("Polygons.txt");


  ASSERT_TRUE(ifs);
  ASSERT_TRUE(ifs.good());

  const int start_num = 1;
  const int end_num = 1000;

  int test_number = start_num;
  while (test_number <= end_num)
  {
    Clipper2Lib::Paths64 subject, subject_open, clip;
    Clipper2Lib::Paths64 solution, solution_open;
    Clipper2Lib::ClipType ct;
    Clipper2Lib::FillRule fr;
    int64_t stored_area, stored_count;

    if (!LoadTestNum(ifs, test_number, 
      subject, subject_open, clip, stored_area, stored_count, ct, fr)) break;

    // check Paths64 solutions
    Clipper2Lib::Clipper64 c;
    c.AddSubject(subject);
    c.AddOpenSubject(subject_open);
    c.AddClip(clip);
    c.Execute(ct, fr, solution, solution_open);

    const int64_t measured_area = static_cast<int64_t>(Area(solution));
    const int64_t measured_count = static_cast<int64_t>(solution.size() + solution_open.size());    
    const int64_t count_diff = stored_count <= 0 ? 0 : std::abs(measured_count - stored_count);
    const int64_t area_diff = stored_area <= 0 ? 0 : std::abs(measured_area - stored_area);
    double area_diff_ratio = (area_diff == 0) ? 0 : std::fabs((double)(area_diff) / measured_area);
    
    // check the polytree variant too
    Clipper2Lib::PolyTree64 solution_polytree;
    Clipper2Lib::Paths64 solution_polytree_open;
    Clipper2Lib::Clipper64 clipper_polytree;
    clipper_polytree.AddSubject(subject);
    clipper_polytree.AddOpenSubject(subject_open);
    clipper_polytree.AddClip(clip);
    clipper_polytree.Execute(ct, fr, solution_polytree, solution_polytree_open);
    
    const int64_t measured_area_polytree = 
      static_cast<int64_t>(solution_polytree.Area());
    const auto solution_polytree_paths = PolyTreeToPaths64(solution_polytree);
    const int64_t measured_count_polytree = 
      static_cast<int64_t>(solution_polytree_paths.size());

    // check polygon counts
    if (stored_count <= 0)
      ; // skip count
    else if (IsInList(test_number, { 120, 121, 130, 138,
      140, 148, 163, 165, 166, 167, 168, 172, 175, 178, 180 }))
      EXPECT_LE(count_diff, 5) << " in test " << test_number;
    else if (IsInList(test_number, { 27, 181 }))
      EXPECT_LE(count_diff, 2) << " in test " << test_number;
    else if (test_number >= 120 && test_number <= 184)
      EXPECT_LE(count_diff, 2) << " in test " << test_number;
    else if (IsInList(test_number, { 23, 45, 87, 102, 111, 113, 191 }))
      EXPECT_LE(count_diff, 1) << " in test " << test_number;
    else
      EXPECT_EQ(count_diff, 0) << " in test " << test_number;

    // check polygon areas
    if (stored_area <= 0)
      ; // skip area
    else if (IsInList(test_number, { 19, 22, 23, 24 }))
      EXPECT_LE((double)area_diff_ratio, 0.5) << " in test " << test_number;
    else if (test_number == 193)
      EXPECT_LE((double)area_diff_ratio, 0.2) << " in test " << test_number;
    else if (test_number == 63)
      EXPECT_LE((double)area_diff_ratio, 0.1) << " in test " << test_number;
    else if (test_number == 16)
      EXPECT_LE((double)area_diff_ratio, 0.075) << " in test " << test_number;
    else if (test_number == 26)
      EXPECT_LE((double)area_diff_ratio, 0.05) << " in test " << test_number;
    else if (IsInList(test_number, { 15, 52, 53, 54, 59, 60, 64, 117, 119, 184 }))
      EXPECT_LE((double)area_diff_ratio, 0.02) << " in test " << test_number;
    else
      EXPECT_LE((double)area_diff_ratio, 0.01) << " in test " << test_number;

    EXPECT_EQ(measured_count, measured_count_polytree) 
      << " in test " << test_number;
    EXPECT_EQ(measured_area, measured_area_polytree) 
      << " in test " << test_number;

    ++test_number;
  }
  //EXPECT_GE(test_number, 188);

  Clipper2Lib::PathsD subjd, clipd, solutiond;
  Clipper2Lib::FillRule frd = Clipper2Lib::FillRule::NonZero;

}