/**
 * @file   StopWatch_test.cc
 * @brief  Test for the StopWatch class
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 3, 2016
 * @see    StopWatch.h
 *
 * This test instantiates a stopwatch and makes some queries.
 *
 * It is not an actual test in that there is no automatic verification of the
 * results (which would be somehow problematic given that times are not
 * completely reproducible).
 *
 */

// our libraries
#include "larcorealg/TestUtils/StopWatch.h"

// C/C++ standard libraries
#include <chrono>
#include <ratio>
#include <thread> // std::this_thread
#include <iostream>


void wait(std::chrono::milliseconds dur) {
  std::cout << " <waiting for " << dur.count() << " ms>" << std::endl;
  std::this_thread::sleep_for(dur);
} // wait()


int main(int /* argc */, char** /* argv */) {

  std::chrono::milliseconds WaitFor { 250 };

  std::cout << "Creating stop watch..." << std::endl;
  testing::StopWatch<std::chrono::duration<double>, std::chrono::system_clock> timer;
  std::cout << " - elapsed time so far: " << timer.elapsed() << " s; partial time: " << timer.partial() << " s" << std::endl;

  wait(WaitFor);
  std::cout << " - elapsed time so far: " << timer.elapsed() << " s (" << timer.elapsed<std::micro>() << " us); partial time: " << timer.partial<std::micro>() << " us" << std::endl;

  std::cout << "Stopping watch" << std::endl;
  timer.stop();
  wait(WaitFor);
  std::cout << " - elapsed time so far: " << timer.elapsed() << " s (" << timer.elapsed<std::chrono::microseconds>() << " us); partial time: " << timer.partial<std::chrono::microseconds>() << " us" << std::endl;

  std::cout << "Resuming watch" << std::endl;
  timer.resume();
  wait(WaitFor);
  std::cout << " - elapsed time so far: " << timer.elapsed() << " s (" << timer.elapsed<std::chrono::seconds>() << " s); partial time: " << timer.partial<std::chrono::seconds>() << " s"  << std::endl;

  std::cout << "Restarting watch" << std::endl;
  timer.restart();
  wait(WaitFor);
  std::cout << " - elapsed time so far: " << timer.elapsed() << " s (" << timer.elapsed<std::chrono::milliseconds>() << " ms); partial time: " << timer.partial<std::chrono::milliseconds>() << " ms"  << std::endl;


  return 0;
} // main()
