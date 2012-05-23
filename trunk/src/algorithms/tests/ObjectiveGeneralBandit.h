// -*- Mode: c++ -*-

/** A multi-objective bandit for vector outcomes */
class ObjectiveGeneralBandit
{
public:
	const int n_actions;
	const int n_outcomes;
	RandomNumberGenerator& rng;
	ObjectiveGeneralBandit(int n_actions_, 
                           int n_outcomes_, 
                           RandomNumberGenerator& rng_)
		: n_actions(n_actions_), n_outcomes(n_outcomes_), rng(rng_)
	{
	}
	virtual ~ObjectiveGeneralBandit()
	{}
	/// get the greedy action
	virtual int getGreedyAction(const Vector& payoff) const = 0;
	/// act
	virtual int Act(const Vector& payoff, const Vector& outcome) = 0;
};



class EpsilonGreedyObjectiveGeneralBandit : public ObjectiveGeneralBandit
{
protected:
	int action;
public:
	Matrix N; ///< probability matrix
	Matrix P; ///< probability matrix
	real epsilon; ///< randomness
	EpsilonGreedyObjectiveGeneralBandit(int n_actions_,
                                        int n_outcomes_,
                                        RandomNumberGenerator& rng_,
                                        real epsilon_)
		: ObjectiveGeneralBandit(n_actions_, n_outcomes_, rng_),
		  action(-1),
		  N(n_actions, n_outcomes),
		  P(n_actions, n_outcomes),
		  epsilon(epsilon_)
	{
		real p = 1.0 / (real) n_outcomes;
		for (int i=0; i<n_actions; ++i) {
			for (int j=0; j<n_outcomes; ++j) {
				N(i, j) = 0.5;
				P(i, j) = p;
			}
		}
		assert(epsilon >= 0.0 && epsilon <= 1.0);
	}

	virtual ~EpsilonGreedyObjectiveGeneralBandit()
	{
	}

	/// get the greedy action
	virtual int getGreedyAction(const Vector& payoff) const
	{
		const Matrix& rP = P;
		return ArgMax(rP * payoff);
	}

	/// act
	virtual int Act(const Vector& payoff, const Vector& outcome)
	{
		// use previous action
		if (action >= 0) {
			Vector p = N.getRow(action) + outcome;
            N.setRow(action, p);
			P.setRow(action, p / p.Sum());
		}

		if (rng.uniform() < epsilon) {
			action = rng.discrete_uniform(n_actions);
		} else {
			action = getGreedyAction(payoff);
		}
		return action;
	}
};


class HoeffdingObjectiveGeneralBandit : public ObjectiveGeneralBandit
{
protected:
	int action;
public:
	Matrix N; ///< incidence matrix
	Matrix P; ///< probability matrix
	real delta; ///< randomness
	int T; ///< number of time-steps
	std::vector<int> plays; ///< number of plays
	HoeffdingObjectiveGeneralBandit(int n_actions_,
                                    int n_outcomes_,
                                    RandomNumberGenerator& rng_)
		: ObjectiveGeneralBandit(n_actions_, n_outcomes_, rng_),
		  action(-1),
		  N(n_actions, n_outcomes),
		  P(n_actions, n_outcomes),
		  delta(1.0),
		  T(0),
		  plays(n_actions)
	{
		real p = 1.0 / (real) n_outcomes;
		for (int i=0; i<n_actions; ++i) {
			plays[i] = 0;
			for (int j=0; j<n_outcomes; ++j) {
				N(i,j) = 0.5;
				P(i,j) = p;
			}
		}
	}

	virtual ~HoeffdingObjectiveGeneralBandit()
	{
	}

    /// Here we have no constraints, so we just add epsilon to all the bits.
	Vector OptimisticTransition(const Vector& p, const Vector& r, real epsilon)
	{		
        int m = p.Size();
        Vector hp = p + epsilon;
        for (int i=0; i<m; ++i) {
            if (hp(i) > 1.0) {
                hp(i) = 1.0;
            }
        }
		return hp;
	}


	/// get the greedy action
	virtual int getGreedyAction(const Vector& payoff) const
	{
		const Matrix& rP = P;
		return ArgMax(rP * payoff);
	}

	/// act
	virtual int Act(const Vector& payoff, const Vector& outcome)
	{
		// use previous action to update probabilities
		if (action >= 0) {
            plays[action]++;
			Vector p = N.getRow(action) + outcome;
            N.setRow(action, p);
			P.setRow(action, p / p.Sum());
		}

		// choose next action
		real u_max = - 1;
		T++;
        delta = 0.1;
		//delta = 1.0 / (real) T;
		for (int i=0; i<n_actions; ++i) {
			real epsilon = HoeffdingBound(n_outcomes, plays[i], delta);
			Vector B = OptimisticTransition(P.getRow(i), payoff, epsilon);
			real u = Product(B, payoff);
			if (u > u_max) {
				u_max = u;
				action = i;
			}
		}

		return action;
	}
};


