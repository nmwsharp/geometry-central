// Copyright (C) 2008 Danil Kirsanov, MIT License
// (Modified to work in geometry-central. Original code can be found here: https://code.google.com/p/geodesic/)

namespace geometrycentral {
namespace surface {

template <class T>
void MemoryAllocator<T>::clear() {
  reset(m_block_size, m_max_number_of_blocks);
}

template <class T>
void MemoryAllocator<T>::reset(unsigned block_size, unsigned max_number_of_blocks) {
  m_block_size = block_size;
  m_max_number_of_blocks = max_number_of_blocks;

  assert(m_block_size > 0);
  assert(m_max_number_of_blocks > 0);

  m_current_position = 0;

  m_storage.reserve(max_number_of_blocks);
  m_storage.resize(1);
  m_storage[0].resize(block_size);

  m_deleted.clear();
  m_deleted.reserve(2 * block_size);
};

template <class T>
typename MemoryAllocator<T>::pointer MemoryAllocator<T>::allocate() {
  pointer result;
  if (m_deleted.empty()) {
    if (m_current_position + 1 >= m_block_size) {
      m_storage.push_back(std::vector<T>());
      m_storage.back().resize(m_block_size);
      m_current_position = 0;
    }
    result = &m_storage.back()[m_current_position];
    ++m_current_position;
  } else {
    result = m_deleted.back();
    m_deleted.pop_back();
  }

  return result;
};

template <class T>
void MemoryAllocator<T>::deallocate(pointer p) {
  if (m_deleted.size() < m_deleted.capacity()) {
    m_deleted.push_back(p);
  }
};

} // namespace surface
} // namespace geometrycentral
