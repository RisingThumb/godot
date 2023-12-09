// Copyright 2014 The Crashpad Authors
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

#include "client/simulate_crash.h"

#include <mach/mach.h>
#include <string.h>
#include <sys/types.h>

#include <iterator>

#include "base/strings/stringprintf.h"
#include "build/build_config.h"
#include "gtest/gtest.h"
#include "test/mac/mach_errors.h"
#include "test/mac/mach_multiprocess.h"
#include "util/mach/exc_server_variants.h"
#include "util/mach/exception_behaviors.h"
#include "util/mach/exception_ports.h"
#include "util/mach/mach_extensions.h"
#include "util/mach/mach_message.h"
#include "util/mach/mach_message_server.h"
#include "util/mach/symbolic_constants_mach.h"
#include "util/misc/implicit_cast.h"

namespace crashpad {
namespace test {
namespace {

class TestSimulateCrashMac final : public MachMultiprocess,
                                   public UniversalMachExcServer::Interface {
 public:
  // Defines which targets the child should set an EXC_CRASH exception handler
  // for.
  enum ExceptionPortsTarget {
    // The child should clear its EXC_CRASH handler for both its task and thread
    // targets. SimulateCrash() will attempt to deliver the exception to the
    // host target, which will fail if not running as root. In any case, the
    // parent should not expect to receive any exception message from the child.
    kExceptionPortsTargetNone = 0,

    // The child will set an EXC_CRASH handler for its task target, and clear it
    // for its thread target. The parent runs an exception server to receive
    // the child’s simulated crash message.
    kExceptionPortsTargetTask,

    // The child will set an EXC_CRASH handler for its thread target, and clear
    // it for its task target. The parent runs an exception server to receive
    // the child’s simulated crash message.
    kExceptionPortsTargetThread,

    // The child sets an EXC_CRASH handler for both its task and thread targets.
    // The parent runs an exception server to receive the message expected to be
    // delivered to the thread target, but returns an error code. The child will
    // then fall back to trying the server registered for the task target,
    // sending a second message to the parent. The server in the parent will
    // handle this one successfully.
    kExceptionPortsTargetBoth,
  };

  TestSimulateCrashMac(ExceptionPortsTarget target,
                       exception_behavior_t behavior,
                       thread_state_flavor_t flavor)
      : MachMultiprocess(),
        UniversalMachExcServer::Interface(),
        target_(target),
        behavior_(behavior),
        flavor_(flavor),
        succeed_(true) {
  }

  TestSimulateCrashMac(const TestSimulateCrashMac&) = delete;
  TestSimulateCrashMac& operator=(const TestSimulateCrashMac&) = delete;

  ~TestSimulateCrashMac() {}

  // UniversalMachExcServer::Interface:
  kern_return_t CatchMachException(exception_behavior_t behavior,
                                   exception_handler_t exception_port,
                                   thread_t thread,
                                   task_t task,
                                   exception_type_t exception,
                                   const mach_exception_data_type_t* code,
                                   mach_msg_type_number_t code_count,
                                   thread_state_flavor_t* flavor,
                                   ConstThreadState old_state,
                                   mach_msg_type_number_t old_state_count,
                                   thread_state_t new_state,
                                   mach_msg_type_number_t* new_state_count,
                                   const mach_msg_trailer_t* trailer,
                                   bool* destroy_complex_request) override {
    *destroy_complex_request = true;

    // Check the entire exception message, because most or all of it was
    // generated by SimulateCrash() instead of the kernel.

    EXPECT_EQ(behavior, behavior_);
    EXPECT_EQ(exception_port, LocalPort());
    if (ExceptionBehaviorHasIdentity(behavior)) {
      EXPECT_NE(thread, THREAD_NULL);
      EXPECT_EQ(task, ChildTask());
    } else {
      EXPECT_EQ(thread, THREAD_NULL);
      EXPECT_EQ(task, TASK_NULL);
    }
    EXPECT_EQ(exception, kMachExceptionSimulated);
    EXPECT_EQ(code_count, 2u);
    if (code_count >= 1) {
      EXPECT_EQ(code[0], 0);
    }
    if (code_count >= 2) {
      EXPECT_EQ(code[1], 0);
    }
    if (!ExceptionBehaviorHasState(behavior)) {
      EXPECT_EQ(*flavor, THREAD_STATE_NONE);
    } else {
      EXPECT_EQ(*flavor, flavor_);
      switch (*flavor) {
#if defined(ARCH_CPU_X86_FAMILY)
        case x86_THREAD_STATE: {
          EXPECT_EQ(old_state_count, x86_THREAD_STATE_COUNT);
          const x86_thread_state* state =
              reinterpret_cast<const x86_thread_state*>(old_state);
          switch (state->tsh.flavor) {
            case x86_THREAD_STATE32:
              EXPECT_EQ(implicit_cast<uint32_t>(state->tsh.count),
                        implicit_cast<uint32_t>(x86_THREAD_STATE32_COUNT));
              break;
            case x86_THREAD_STATE64:
              EXPECT_EQ(implicit_cast<uint32_t>(state->tsh.count),
                        implicit_cast<uint32_t>(x86_THREAD_STATE64_COUNT));
              break;
            default:
              ADD_FAILURE() << "unexpected tsh.flavor " << state->tsh.flavor;
              break;
          }
          break;
        }
        case x86_FLOAT_STATE: {
          EXPECT_EQ(old_state_count, x86_FLOAT_STATE_COUNT);
          const x86_float_state* state =
              reinterpret_cast<const x86_float_state*>(old_state);
          switch (state->fsh.flavor) {
            case x86_FLOAT_STATE32:
              EXPECT_EQ(implicit_cast<uint32_t>(state->fsh.count),
                        implicit_cast<uint32_t>(x86_FLOAT_STATE32_COUNT));
              break;
            case x86_FLOAT_STATE64:
              EXPECT_EQ(implicit_cast<uint32_t>(state->fsh.count),
                        implicit_cast<uint32_t>(x86_FLOAT_STATE64_COUNT));
              break;
            default:
              ADD_FAILURE() << "unexpected fsh.flavor " << state->fsh.flavor;
              break;
          }
          break;
        }
        case x86_DEBUG_STATE: {
          EXPECT_EQ(old_state_count, x86_DEBUG_STATE_COUNT);
          const x86_debug_state* state =
              reinterpret_cast<const x86_debug_state*>(old_state);
          switch (state->dsh.flavor) {
            case x86_DEBUG_STATE32:
              EXPECT_EQ(implicit_cast<uint32_t>(state->dsh.count),
                        implicit_cast<uint32_t>(x86_DEBUG_STATE32_COUNT));
              break;
            case x86_DEBUG_STATE64:
              EXPECT_EQ(implicit_cast<uint32_t>(state->dsh.count),
                        implicit_cast<uint32_t>(x86_DEBUG_STATE64_COUNT));
              break;
            default:
              ADD_FAILURE() << "unexpected dsh.flavor " << state->dsh.flavor;
              break;
          }
          break;
        }
        case x86_THREAD_STATE32:
          EXPECT_EQ(old_state_count, x86_THREAD_STATE32_COUNT);
          break;
        case x86_FLOAT_STATE32:
          EXPECT_EQ(old_state_count, x86_FLOAT_STATE32_COUNT);
          break;
        case x86_DEBUG_STATE32:
          EXPECT_EQ(old_state_count, x86_DEBUG_STATE32_COUNT);
          break;
        case x86_THREAD_STATE64:
          EXPECT_EQ(old_state_count, x86_THREAD_STATE64_COUNT);
          break;
        case x86_FLOAT_STATE64:
          EXPECT_EQ(old_state_count, x86_FLOAT_STATE64_COUNT);
          break;
        case x86_DEBUG_STATE64:
          EXPECT_EQ(old_state_count, x86_DEBUG_STATE64_COUNT);
          break;
#elif defined(ARCH_CPU_ARM64)
        case ARM_UNIFIED_THREAD_STATE: {
          EXPECT_EQ(old_state_count, ARM_UNIFIED_THREAD_STATE_COUNT);
          const arm_unified_thread_state* state =
              reinterpret_cast<const arm_unified_thread_state*>(old_state);
          EXPECT_EQ(state->ash.flavor,
                    implicit_cast<uint32_t>(ARM_THREAD_STATE64));
          if (state->ash.flavor == ARM_THREAD_STATE64) {
            EXPECT_EQ(state->ash.count,
                      implicit_cast<uint32_t>(ARM_THREAD_STATE64_COUNT));
          }
          break;
        }
        case ARM_THREAD_STATE64:
          EXPECT_EQ(old_state_count, ARM_THREAD_STATE64_COUNT);
          break;
        case ARM_NEON_STATE64:
          EXPECT_EQ(old_state_count, ARM_NEON_STATE64_COUNT);
          break;
        case ARM_DEBUG_STATE64:
          EXPECT_EQ(old_state_count, ARM_DEBUG_STATE64_COUNT);
          break;
#else
#error Port to your CPU architecture
#endif
        default:
          ADD_FAILURE() << "unexpected flavor " << *flavor;
          break;
      }

      // Attempt to set a garbage thread state, which would cause the child to
      // crash inside SimulateCrash() if it actually succeeded. This tests that
      // SimulateCrash() ignores new_state instead of attempting to set the
      // state as the kernel would do. This operates in conjunction with the
      // |true| argument to ExcServerSuccessfulReturnValue() below.
      *new_state_count = old_state_count;
      size_t new_state_size = sizeof(natural_t) * old_state_count;
      memset(new_state, 0xa5, new_state_size);
    }

    if (!succeed_) {
      // The client has registered EXC_CRASH handlers for both its thread and
      // task targets, and sent a simulated exception message to its
      // thread-level EXC_CRASH handler. To test that it will fall back to
      // trying the task-level EXC_CRASH handler, return a failure code, which
      // should cause SimulateCrash() to try the next target.
      EXPECT_EQ(target_, kExceptionPortsTargetBoth);
      return KERN_ABORTED;
    }

    ExcServerCopyState(
        behavior, old_state, old_state_count, new_state, new_state_count);

    return ExcServerSuccessfulReturnValue(exception, behavior, true);
  }

 private:
  // MachMultiprocess:

  void MachMultiprocessParent() override {
    if (target_ == kExceptionPortsTargetNone) {
      // The child does not have any EXC_CRASH handlers registered for its
      // thread or task targets, so no exception message is expected to be
      // generated. Don’t run the server at all.
      return;
    }

    UniversalMachExcServer universal_mach_exc_server(this);

    mach_msg_return_t mr;
    if (target_ == kExceptionPortsTargetBoth) {
      // The client has registered EXC_CRASH handlers for both its thread and
      // task targets. Run a server that will return a failure code when the
      // exception message is sent to the thread target, which will cause the
      // client to fall back to the task target and send another message.
      succeed_ = false;
      mr = MachMessageServer::Run(&universal_mach_exc_server,
                                  LocalPort(),
                                  MACH_MSG_OPTION_NONE,
                                  MachMessageServer::kOneShot,
                                  MachMessageServer::kReceiveLargeError,
                                  kMachMessageTimeoutWaitIndefinitely);
      EXPECT_EQ(mr, MACH_MSG_SUCCESS)
          << MachErrorMessage(mr, "MachMessageServer::Run");
    }

    succeed_ = true;
    mr = MachMessageServer::Run(&universal_mach_exc_server,
                                LocalPort(),
                                MACH_MSG_OPTION_NONE,
                                MachMessageServer::kOneShot,
                                MachMessageServer::kReceiveLargeError,
                                kMachMessageTimeoutWaitIndefinitely);
    EXPECT_EQ(mr, MACH_MSG_SUCCESS)
        << MachErrorMessage(mr, "MachMessageServer::Run");
  }

  void MachMultiprocessChild() override {
    bool task_valid = target_ == kExceptionPortsTargetTask ||
                      target_ == kExceptionPortsTargetBoth;
    ExceptionPorts task_exception_ports(ExceptionPorts::kTargetTypeTask,
                                        TASK_NULL);
    ASSERT_TRUE(task_exception_ports.SetExceptionPort(
        EXC_MASK_CRASH,
        task_valid ? RemotePort() : MACH_PORT_NULL,
        behavior_,
        flavor_));

    bool thread_valid = target_ == kExceptionPortsTargetThread ||
                        target_ == kExceptionPortsTargetBoth;
    ExceptionPorts thread_exception_ports(ExceptionPorts::kTargetTypeThread,
                                          THREAD_NULL);
    ASSERT_TRUE(thread_exception_ports.SetExceptionPort(
        EXC_MASK_CRASH,
        thread_valid ? RemotePort() : MACH_PORT_NULL,
        behavior_,
        flavor_));

    CRASHPAD_SIMULATE_CRASH();
  }

  ExceptionPortsTarget target_;
  exception_behavior_t behavior_;
  thread_state_flavor_t flavor_;
  bool succeed_;
};

TEST(SimulateCrash, SimulateCrash) {
  static constexpr TestSimulateCrashMac::ExceptionPortsTarget kTargets[] = {
      TestSimulateCrashMac::kExceptionPortsTargetNone,
      TestSimulateCrashMac::kExceptionPortsTargetTask,
      TestSimulateCrashMac::kExceptionPortsTargetThread,
      TestSimulateCrashMac::kExceptionPortsTargetBoth,
  };

  static constexpr exception_behavior_t kBehaviors[] = {
      EXCEPTION_DEFAULT,
      EXCEPTION_STATE,
      EXCEPTION_STATE_IDENTITY,
      EXCEPTION_DEFAULT | kMachExceptionCodes,
      EXCEPTION_STATE | kMachExceptionCodes,
      EXCEPTION_STATE_IDENTITY | kMachExceptionCodes,
  };

  static constexpr thread_state_flavor_t kFlavors[] = {
#if defined(ARCH_CPU_X86_FAMILY)
    x86_THREAD_STATE,
    x86_FLOAT_STATE,
    x86_DEBUG_STATE,
#if defined(ARCH_CPU_X86)
    x86_THREAD_STATE32,
    x86_FLOAT_STATE32,
    x86_DEBUG_STATE32,
#elif defined(ARCH_CPU_X86_64)
    x86_THREAD_STATE64,
    x86_FLOAT_STATE64,
    x86_DEBUG_STATE64,
#endif
#elif defined(ARCH_CPU_ARM64)
    ARM_UNIFIED_THREAD_STATE,
    ARM_THREAD_STATE64,
    ARM_NEON_STATE64,
    ARM_DEBUG_STATE64,
#else
#error Port to your CPU architecture
#endif
  };

  for (size_t target_index = 0; target_index < std::size(kTargets);
       ++target_index) {
    TestSimulateCrashMac::ExceptionPortsTarget target = kTargets[target_index];
    SCOPED_TRACE(base::StringPrintf(
        "target_index %zu, target %d", target_index, target));

    for (size_t behavior_index = 0; behavior_index < std::size(kBehaviors);
         ++behavior_index) {
      exception_behavior_t behavior = kBehaviors[behavior_index];
      SCOPED_TRACE(base::StringPrintf(
          "behavior_index %zu, behavior %s",
          behavior_index,
          ExceptionBehaviorToString(behavior, kUseFullName | kUnknownIsNumeric)
              .c_str()));

      if (!ExceptionBehaviorHasState(behavior)) {
        TestSimulateCrashMac test_simulate_crash_mac(
            target, behavior, THREAD_STATE_NONE);
        test_simulate_crash_mac.Run();
      } else {
        for (size_t flavor_index = 0; flavor_index < std::size(kFlavors);
             ++flavor_index) {
          thread_state_flavor_t flavor = kFlavors[flavor_index];
          SCOPED_TRACE(base::StringPrintf(
              "flavor_index %zu, flavor %s",
              flavor_index,
              ThreadStateFlavorToString(
                  flavor, kUseFullName | kUnknownIsNumeric).c_str()));

          TestSimulateCrashMac test_simulate_crash_mac(
              target, behavior, flavor);
          test_simulate_crash_mac.Run();
        }
      }
    }
  }
}

}  // namespace
}  // namespace test
}  // namespace crashpad
