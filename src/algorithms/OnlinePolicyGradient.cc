#include "OnlinePolicyGradient.h"

PolicyGradientActorCritic::PolicyGradientActorCritic(int n_states_, int n_actions_, real gamma_, real step_size_) :
	n_states(n_states_),
	n_actions(n_actions_),
	gamma(gamma_),
	step_size(step_size_),
	critic(n_states, n_actions, gamma, 0.0, step_size),
	policy(n_states, n_actions),
	params(n_states, n_actions),
	Q(n_states, n_actions)
{
	Reset();
}

real PolicyGradientActorCritic::Observe(real reward, const int& next_state, const int& next_action)
{
	real d = 0;
	if (state >= 0) {
		d = GradientUpdate(state, action); // not sure how much difference it makes to update the next state-action pair instead
		//printf("%d %d %d %d%f\n", state, action, d);
		Q(state, action) += step_size * (reward + gamma * Q(next_state, next_action) - Q(state, action));
	}
	critic.Observe(reward, next_state, next_action);
	UpdatePolicy();
	state = next_state;
	action = next_action;
	return d;
}

int PolicyGradientActorCritic::Act(real reward, const int& next_state)
{
	policy.Observe(reward, next_state);
	int next_action = policy.SelectAction();
	Observe(reward, next_state, next_action);
	return next_action;
}

real PolicyGradientActorCritic::GradientUpdate(int s, int a)
{
	//real U = critic.getValue(s);
	//real U = critic.getValue(s);
	real U = Q(s,a);
	printf("s:%d, a:%d: Q_w:%f Q_S:%f\n", s, a, U, critic.getValue(s,a));
	real delta = 0;
	for (int j=0; j<n_actions; ++j) {
		real p_a = policy.getActionProbability(s, j);
		real d_sj =  0;
		if (j==a) {
			d_sj = (1.0 - p_a) * U;
		} else {
			d_sj = -p_a * U;
		}
		params(s,j) += step_size * d_sj;
		delta += fabs(d_sj);
	}
	return delta;
}

void PolicyGradientActorCritic::UpdatePolicy()
{
	for (int s=0; s<n_states; ++s) {
		Vector eW = exp(params.getRow(s));
		eW /= eW.Sum();
		Vector* pS = policy.getActionProbabilitiesPtr(s);
		(*pS) = eW;
	}
}

// --- PGAC with features --- //

PolicyGradientActorCriticPhi::PolicyGradientActorCriticPhi(BasisSet<Vector, int>& basis_,
														   int n_states_,
														   int n_actions_,
														   real gamma_,
														   real step_size_) :
	basis(basis_),
	n_states(n_states_),
	n_actions(n_actions_),
	gamma(gamma_),
	step_size(step_size_),
	critic(n_states, n_actions, basis, gamma),
	policy(n_states, n_actions, basis), // bug?
	params(basis.size())

{
	Reset();
}

real PolicyGradientActorCriticPhi::Observe(real reward, const Vector& next_state, const int& next_action)
{
	basis.Observe(action, reward, next_state);
	real d = 0;
	if (valid_state) {
		d = GradientUpdate(state, action);
	}
	//critic.Observe(reward, next_state, next_action);
	UpdatePolicy();
	state = next_state;
	action = next_action;
	valid_state = true;
	return d;
}

int PolicyGradientActorCriticPhi::Act(real reward, const Vector& next_state)
{
	basis.Observe(action, reward, next_state);
	Vector features = basis.F();
	int next_action = policy.SelectAction(next_state);
	Observe(reward, next_state, next_action);
	return next_action;
}

real PolicyGradientActorCriticPhi::GradientUpdate(const Vector& s, int a)
{
	basis.Observe(s);
	Vector phi = basis.F();
	// copy the state-features into a larger state-action feature vector
	Vector features(phi.Size()*n_actions);
	int k = a * phi.Size();
	for (int i=0; i<phi.Size(); ++i) {
		features(k) = phi(i);
	}
	real U = critic.getValue(s);	
	policy.GradientUpdate(s, a, U);
	return 0;
}


void PolicyGradientActorCriticPhi::UpdatePolicy()
{

}
