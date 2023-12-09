// Copyright 2020 The Crashpad Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#import <UIKit/UIKit.h>
#import <XCTest/XCTest.h>

#include "base/check.h"
#import "test/ios/cptest_google_test_runner_delegate.h"

#if !defined(__has_feature) || !__has_feature(objc_arc)
#error "This file requires ARC support."
#endif

@interface CPTestGoogleTestRunner : XCTestCase
@end

@implementation CPTestGoogleTestRunner

- (void)testRunGoogleTests {
  id appDelegate = UIApplication.sharedApplication.delegate;
  DCHECK([appDelegate
      conformsToProtocol:@protocol(CPTestGoogleTestRunnerDelegate)]);

  id<CPTestGoogleTestRunnerDelegate> runnerDelegate =
      static_cast<id<CPTestGoogleTestRunnerDelegate>>(appDelegate);
  DCHECK(runnerDelegate.supportsRunningGoogleTestsWithXCTest);
  XCTAssertTrue([runnerDelegate runGoogleTests] == 0);
}

@end
