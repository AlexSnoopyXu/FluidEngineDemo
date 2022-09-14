#include "AnimationBase.h"

void AnimationBase::update(const Frame& frame)
{
	//Some pre-processing here

	onUpdate(frame);

	//Some post-processing here
}

void PhysicsAnimation::onUpdate(const Frame& frame) 
{
	if (frame.index > _currentFrame.index)
	{
		unsigned int numOfFrames = frame.index - _currentFrame.index;
		for (size_t i = 0; i < numOfFrames; i++)
		{
			advanceTimeStep(frame.timeIntervalInSeconds);
		}

		_currentFrame = frame;
	}
}

void PhysicsAnimation::advanceTimeStep(double timeIntervalInSeconds)
{
	onAdvanceTimStep(timeIntervalInSeconds);
}