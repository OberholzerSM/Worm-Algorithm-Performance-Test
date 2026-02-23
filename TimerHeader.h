#pragma once
#include <chrono>

class Timer
{
public:
	Timer()
	{
		start_time = clock::now();
	}

	void reset()
	{
		start_time = clock::now();
		pause_time = start_time;
	}

	void pause()
	{
		if(!paused)
			pause_time = clock::now();
		paused = true;
	}

	void resume()
	{
		if (paused)
			start_time += clock::now() - pause_time;
		paused = false;
	}

	double getTime() const
	{
		if (!paused)
			return std::chrono::duration_cast<double_second>(clock::now() - start_time).count();
		else
			return std::chrono::duration_cast<double_second>(pause_time - start_time).count();
	}

private:
	using clock = std::chrono::steady_clock;								//The clock-type to use.
	using time_point = std::chrono::time_point<clock>;						//Defined a single point in time using the clock.
	using double_second = std::chrono::duration<double, std::ratio<1> >;	//Defined the Type seconds using double.
	time_point start_time, pause_time;										//When the clock starts counting and when the clock got paused.
	bool paused{false};
};