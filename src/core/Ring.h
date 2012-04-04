/* -*- Mode: c++;  -*- */
// copyright (c) 2009 by Christos Dimitrakakis <christos.dimitrakakis@gmail.com>
/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#ifndef RING_H
#define RING_H

#include <vector>
#include <iostream>
#include <cassert>
#include "SmartAssert.h"

// forward declaration
//template <typename C>
//class Ring;

/** A ring buffer.

    A ring buffer is useful for many things.
    
 */
template <typename C>
class Ring
{
public:
    class iterator
    {
    protected:
        const Ring<C>* ring;
        int index;
        int length;
    public:
        iterator()
        {
        }
        iterator(Ring<C>* ring_)
            : ring(ring_),
              index(ring->pos),
              length(ring->data.size())
        {
        }
        void operator--()
        {
            index++;
            while (length && index >= length) {
                index -= length;
            }
        }
        void operator++()
        {
            index--;
           while (length && index < 0) {
                index += length;
            }
        }
		void SetIndex(int i)
		{
			assert(length == 0 || (i >=0 && i < length));
			index = i;
		}
		int GetIndex()
		{
			return index;
		}
        int operator*()
        {
            assert(index >=0 && index < length);
            return ring->data[index];
        }
        bool operator==(const iterator& rhs) const
        {
            return (ring == rhs.ring && index == rhs.index);
        }
        bool operator!=(const iterator& rhs) const
        {
            return (ring != rhs.ring || index != rhs.index);
        }
           
    };
protected:
    C null_element;
    int max_size; ///< the maximum size of the array
    std::vector<C> data; ///< the data
    int pos; ///< the current position of the last write in the ring buffer. 
    ulong T;
    iterator end_iterator;
public:
    /// Construct a buffer of a given size
    Ring(int max_size_) : null_element(0),
                          max_size(max_size_),
                          data(max_size + 1),
                          pos(0),
                          T(0),
                          end_iterator(this)
    {
        clear();
    }


	/// Clear data and reset position and end iterator.
    void clear()
    {
        for (uint i=0; i<data.size(); ++i) {
            data[i] = 0;
        }
        pos = max_size;
        end_iterator.SetIndex(max_size);
    }

    /// Construct a buffer of a given size
    Ring() : null_element(0), max_size(0), pos(0), T(0)
    {
    }

    /// Destructor.
    ~Ring()
    {
    }

	/// Get an iterator to the current position
    iterator begin()
    {
        return (iterator (this));
    }

	/// Get an iterator past the final position
    iterator end()
    {
        return end_iterator;
    }

	/// Change the size of the ring. This clears everything.
    void resize(size_t new_size)
    {
        max_size = new_size;
        data.resize(max_size + 1);
        clear();
    }

    /** Get the unique ID of a state.
        
        \note Only use with care. It can fail due to int overflow.

        A better implementation may exist with hashes.
    */
    long long get_id(int n_states)
    {
        if (max_size ==0) {
            return 0;
        }
        int n = 1;
        long long id = 0;
        for (iterator it = begin(); it != end_iterator; ++it, n*=n_states) {
            id += n * (*it);
        }
        return id;
    }

    /** Push back a datum.
        
        Whenever new data is available, it is inserted in the current
        position. The position pointer wraps around.

        @param x The data to be pushed back.
     */
    void push_back(C x)
    {
        if (max_size) {
            pos++;
            if (pos >= (int) data.size()) {
                pos -= data.size();
            }
            data[pos] = x;

			//std::cout << "Moving front iterator to " << pos << std::endl;
			// only change the end iterator if we've hit it
			if (end_iterator.GetIndex() == pos) {
				int end_pos = pos + 1;
				if (end_pos >= (int) data.size()) {
					end_pos -= data.size();
				}
				end_iterator.SetIndex(end_pos);
				//std::cout << "Moving end iterator to " << end_pos << std::endl;
			}
        }
        T++;
        //printf ("push: %d %d %d/%d\n", pos, end_iterator.GetIndex(), T, data.size());
    }

    /** Pop a datum.
        
        Remove a datum from the head.
     */
    void pop()
    {
        if (end_iterator.GetIndex() == pos) {
            return;
        }
        if (max_size) {
            pos--;
            if (pos < 0) {
                pos += data.size();
            }
            if (end_iterator.GetIndex() == pos) {
				int end_pos = pos - 1;
				if (end_pos <0) {
					end_pos += data.size();
				}
				end_iterator.SetIndex(end_pos);
			}
        }
        T--;

        //printf ("pop: %d %d %d\n", pos, end_iterator.GetIndex(), T);
    }

    
    /** Access data at offset from the current position.
        
        If we are currently at time T, then we access the data
        inserted at time T - offset.
        
        @param offset This must be non-negative and smaller than \c max_size.

        @return A reference to the data.
     */
    C& operator[](int offset)
    {
        assert(offset >= 0 && offset < max_size);
        int x = (pos - offset);
		/* Alternative implementations

		  1:  x - y*(x/y)
		  2: n = x % abs(y); if (n<0) n += abs(y);
          3: ((x % y) + y) & y
		 */
		x %= max_size;
		if (x < 0) 
		  x += max_size;

        //        std::cout << x << std::endl;
        SMART_ASSERT(x >= 0 && x < max_size)(x)(max_size);

        return data[x];
    }
    
    /// Return the oldest element in the ring buffer
    C& front()
    {
        if (max_size) {
            return data[(pos + 1)%max_size];
        } else { 
            return null_element;
        }
    }

    /// Return the last element.
    C& back()
    {
        if (max_size) {
            return data[pos];
        } else { 
            return null_element;
        }
    }        
    ulong size()
    {
		if (T < (ulong) max_size) {
			return T;
		}
		return max_size;
    }
};

#endif
