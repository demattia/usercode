#ifndef TestObjectWrappers_h
#define TestObjectWrappers_h

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include <iostream>
#include <boost/any.hpp>

/**
 * boost::any is used to pass the objects contained in classes derived from BaseObjectWrapper. To use the method
 * we need to define it in the base classe, but the return type depends on the template parameter. Thus we return
 * a boost::any object containing the actual object. The operator== willing to access the object (and knowing its
 * type) can do the conversion back with boost::any_cast<type>(anyObject).
 */

using namespace std;

namespace {
  typedef edm::ParameterSet PSet;
  typedef std::vector<PSet> VPSet;
  typedef VPSet::const_iterator iter_t;
}

/// Empty base class only used to have a polymorphic pointer
class BaseObjectWrapper {
public:
  virtual bool operator==( const BaseObjectWrapper & other ) = 0;
  virtual bool operator!=( const BaseObjectWrapper & other ) { return !(*this == other); }
  virtual BaseObjectWrapper * eval(const PSet & parameters) { return this; }
  virtual boost::any object() const = 0;
};

/**
 * This class will hold the tested object
 * The object must be assignable
 */
template <class T>
class ObjectWrapper : public BaseObjectWrapper
{
public:
  /// Not initializing empty constructor used by derived classes filling the object_ on their own.
  ObjectWrapper() {}
  /// Build from the object
  ObjectWrapper(const T & test) : object_(test) {}
  virtual ~ObjectWrapper() {}
  virtual bool operator==( const BaseObjectWrapper & other )
  {
    cout << "Warning: using the default operator== which always returns false. It must be overridden in derived classes" << endl;
    return false;
  }
//   {
//     return object_ == boost::any_cast<T>(other.object());
//   }
  boost::any object() const { return boost::any(object_); }
protected:
  T object_;
};

/**
 * Define wrappers for basic types
 */

/**
 * It compares only to bools. </br>
 * It can be used for self tests when the operator== of the object does everything internally
 * so that the success of failure of the test can be specified by this bool.
 */
class BoolWrapper : public ObjectWrapper<bool>
{
public:
  BoolWrapper(const bool b) { object_ = b; }
  virtual ~BoolWrapper() {}
//   virtual bool operator==( const BaseObjectWrapper & other )
//   {
//     // This way we use the operator== defined in the other object
//     BaseObjectWrapper * thisPointer = dynamic_cast<BaseObjectWrapper*>(this);
//     if(other == *thisPointer) return true;
//     // if( object_ == boost::any_cast<bool>(other.object()) ) return true;
//     return false;
//   }
};

#endif
