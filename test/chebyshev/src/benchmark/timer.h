
///
/// @file timer.h A timer class to measure elapsed time in milliseconds.
///

#ifndef CHEBYSHEV_TIMER_H
#define CHEBYSHEV_TIMER_H

#include <chrono>


namespace chebyshev {

	namespace benchmark {

		/// @class timer
		/// Timer class to measure elapsed time in milliseconds.
		class timer {
			private:
				std::chrono::time_point<std::chrono::high_resolution_clock> s;

			public:

				/// Constructs the timer storing the current time.
				timer() {
					start();
				}


				/// Start the timer.
				void start() {
					s = std::chrono::high_resolution_clock::now();
				}


				/// Returns the elapsed time since construction or
				/// start of the timer in milliseconds.
				inline long double get() const {

					auto start = std::chrono::time_point_cast<std::chrono::milliseconds>(s)
									.time_since_epoch();

					auto end = std::chrono::time_point_cast<std::chrono::milliseconds>(
								std::chrono::high_resolution_clock::now())
									.time_since_epoch();

					return (long double) (end - start).count();
				}


				/// Returns the elapsed time since construction or
				/// start of the timer in milliseconds.
				/// @see get
				inline long double operator()() {
					return get();
				}
			
		};

	}
}

#endif
