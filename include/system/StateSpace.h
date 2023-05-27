#pragma once
#include <ControlOutputSystem.h>

namespace ct
{
    namespace core
    {
        template <size_t STATE_DIM, size_t CONTROL_DIM, size_t OUTPUT_DIM>
        class StateSpace : public ControlOutputSystem<STATE_DIM, CONTROL_DIM, OUTPUT_DIM>
        {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            using Base = ControlOutputSystem<STATE_DIM, CONTROL_DIM, OUTPUT_DIM>;
            MAKE_TIME_STATE_CONTROL_OUTPUT_TYPE_FROM_BASE

            StateSpace(const state_t &initState = state_t()) : Base(initState)
            {
                A_.setZero();
                B_.setZero();
                C_.setZero();
            }

            // constructor
            StateSpace(const StateMatrix<STATE_DIM> &A,
                       const StateControlMatrix<STATE_DIM, CONTROL_DIM> &B,
                       const OutputStateMatrix<OUTPUT_DIM, STATE_DIM> &C,
                       const state_t &initState = state_t()) : Base(initState), A_(A), B_(B), C_(C) {}
            // copy constructor
            StateSpace(const StateSpace &other) : Base(other), A_(other.A_), B_(other.B_), C_(other.C_) {}

            StateSpace *clone() const override
            {
                return new StateSpace(*this); // calls copy constructor
            }
            // destructor
            ~StateSpace() = default;
            // The system dynamics. We override this method which gets called by e.g. the Integrator
            void computeControlledDynamics(const state_t &state,
                                           const time_t &t,
                                           const control_t &control,
                                           state_t &derivative) override
            {
                derivative = A_ * state + B_ * control;
            }

            output_t computeOutput(const state_t &state, const time_t &t = 0.0)
            {
                return C_ * state;
            }

            void setSystemMatrices(const StateMatrix<STATE_DIM> &A,
                                   const StateControlMatrix<STATE_DIM, CONTROL_DIM> &B,
                                   const OutputStateMatrix<OUTPUT_DIM, STATE_DIM> &C)
            {
                A_ = A;
                B_ = B;
                C_ = C;
            }

            void getSystemMatrices(Eigen::Matrix<double, STATE_DIM, STATE_DIM> &A,
                                   Eigen::Matrix<double, STATE_DIM, CONTROL_DIM> &B,
                                   Eigen::Matrix<double, OUTPUT_DIM, STATE_DIM> &C)
            {
                A = A_;
                B = B_;
                C = C_;
            }

        protected:
            StateMatrix<STATE_DIM> A_;
            StateControlMatrix<STATE_DIM, CONTROL_DIM> B_;
            OutputStateMatrix<OUTPUT_DIM, STATE_DIM> C_;
        };

    }
}