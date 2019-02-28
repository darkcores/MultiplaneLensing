#pragma once

#include "composite.h"
#include "image.h"

#include <map>
#include <utility>
#include <memory>

class Multiplane {
  private:
	// apparently std map is already sorted on keys, which is probably useful 
	std::map<double, std::shared_ptr<CompositeLens>> lenses;

  public:
    Multiplane();
};
