#include "particle.h"
#include "parameters.h"
#include "likelihood_estimator_output.h"

Particle::Particle()
{
  this->likelihood_estimator_outputs.resize(0);
}

Particle::Particle(const Parameters &parameters_in)
{
  this->parameters = parameters_in;
  this->likelihood_estimator_outputs.resize(0);
}

Particle::Particle(const Parameters &parameters_in,
                   const std::vector<LikelihoodEstimatorOutput*> &outputs_in)
{
  this->parameters = parameters_in;
  this->likelihood_estimator_outputs = outputs_in;
}

Particle::~Particle()
{
  for (std::vector<LikelihoodEstimatorOutput*>::iterator i=this->likelihood_estimator_outputs.begin();
       i!=this->likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      delete *i;
  }
}

//Copy constructor for the Particle class.
Particle::Particle(const Particle &another)
{
  this->make_copy(another);
}

void Particle::operator=(const Particle &another)
{
  if(this == &another)
    return;

  this->make_copy(another);
}

void Particle::make_copy(const Particle &another)
{
  this->parameters = another.parameters;
  //this->model_and_algorithm = another.model_and_algorithm;

  this->likelihood_estimator_outputs.resize(0);
  this->likelihood_estimator_outputs.reserve(another.likelihood_estimator_outputs.size());
  for (std::vector<LikelihoodEstimatorOutput*>::const_iterator i=another.likelihood_estimator_outputs.begin();
       i!=another.likelihood_estimator_outputs.end();
       ++i)
  {
    if (*i!=NULL)
      this->likelihood_estimator_outputs.push_back((*i)->duplicate());
  }
}

std::ostream& operator<<(std::ostream& os, const Particle &p)
{
  os << p.parameters << std::endl;

  os << "{" << std::endl;

  std::vector<LikelihoodEstimatorOutput*>::const_iterator it;

  for (it=p.likelihood_estimator_outputs.begin();it!=p.likelihood_estimator_outputs.end();++it)
  {
    if (it==p.likelihood_estimator_outputs.begin())
      os << *it << std::endl;
    else
      os << "," << *it << std::endl;
  }

  os << "}";

  return os;
}
