#pragma once
#include <ct/core/core.h> // as usual, include CT

#define MAKE_TIME_STATE_CONTROL_OUTPUT_TYPE_FROM_BASE \
    using time_t = typename Base::time_t;             \
    using state_t = typename Base::state_t;           \
    using control_t = typename Base::control_t;       \
    using output_t = typename Base::output_t;

namespace ct
{
    namespace core
    {

        template <size_t STATE_DIM, size_t CONTROL_DIM, size_t OUTPUT_DIM, typename SCALAR = double>
        class ControlOutputSystem : public ControlledSystem<STATE_DIM, CONTROL_DIM>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using Base = ControlledSystem<STATE_DIM, CONTROL_DIM>;
            using time_t = typename Base::time_t;
            typedef StateVector<STATE_DIM> state_t;
            typedef ControlVector<CONTROL_DIM> control_t;
            typedef OutputVector<OUTPUT_DIM> output_t;

            // constructor
            ControlOutputSystem(const state_t &initState) : Base(), initState_(initState) {}
            // copy constructor
            ControlOutputSystem(const ControlOutputSystem &other) : Base(other), initState_(other.initState_) {}

            virtual ControlOutputSystem *clone() const override = 0;
            // destructor
            ~ControlOutputSystem() = default;

            void setInitState(const state_t &initState) { initState_ = initState; }

            virtual void computeControlledDynamics(const state_t &state,
                                                   const time_t &t,
                                                   const control_t &control,
                                                   state_t &derivative) override = 0;

            virtual output_t computeOutput(const state_t &state, const time_t &t) = 0;

            state_t getInitState() const
            {
                return initState_;
            }

        protected:
            state_t initState_;
        };

    }
}