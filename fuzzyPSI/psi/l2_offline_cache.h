#pragma once

#include "mImt.h"
#include "mPeqt/peqt.h"

#include <filesystem>
#include <fstream>
#include <cstring>
#include <string>
#include <type_traits>

namespace CmpFuzzyPSI::l2cache {

struct Key {
  u64 senderSize;
  u64 receiverSize;
  u64 ssp;
  u64 dim;
  u64 metric;
  u64 delta;
};

constexpr u64 kMagic = 0x4c324f4646434143ULL; // "L2OFFCAC"
constexpr u64 kVersion = 2;

template <typename T>
void writeValue(std::ofstream &out, const T &value) {
  static_assert(std::is_trivially_copyable_v<T>);
  out.write(reinterpret_cast<const char *>(&value), sizeof(T));
  if (!out)
    throw std::runtime_error("failed while writing L2 offline cache");
}

template <typename T>
void readValue(std::ifstream &in, T &value) {
  static_assert(std::is_trivially_copyable_v<T>);
  in.read(reinterpret_cast<char *>(&value), sizeof(T));
  if (!in)
    throw std::runtime_error("failed while reading L2 offline cache");
}

template <typename T>
void writeVector(std::ofstream &out, const std::vector<T> &values) {
  static_assert(std::is_trivially_copyable_v<T>);
  const u64 size = values.size();
  writeValue(out, size);
  if (size)
    out.write(reinterpret_cast<const char *>(values.data()), size * sizeof(T));
  if (!out)
    throw std::runtime_error("failed while writing L2 offline cache vector");
}

template <typename T>
void readVector(std::ifstream &in, std::vector<T> &values) {
  static_assert(std::is_trivially_copyable_v<T>);
  u64 size;
  readValue(in, size);
  values.resize(size);
  if (size)
    in.read(reinterpret_cast<char *>(values.data()), size * sizeof(T));
  if (!in)
    throw std::runtime_error("failed while reading L2 offline cache vector");
}

inline void writeBitVector(std::ofstream &out, const BitVector &bits) {
  const u64 size = bits.size();
  writeValue(out, size);
  if (size)
    out.write(reinterpret_cast<const char *>(bits.data()), bits.sizeBytes());
  if (!out)
    throw std::runtime_error("failed while writing L2 offline cache bit vector");
}

inline void readBitVector(std::ifstream &in, BitVector &bits) {
  u64 size;
  readValue(in, size);
  bits.resize(size);
  if (size)
    in.read(reinterpret_cast<char *>(bits.data()), bits.sizeBytes());
  if (!in)
    throw std::runtime_error("failed while reading L2 offline cache bit vector");
}

inline void writeNestedBytes(std::ofstream &out,
                             const std::vector<std::vector<u8>> &values) {
  const u64 size = values.size();
  writeValue(out, size);
  for (const auto &value : values)
    writeVector(out, value);
}

inline void readNestedBytes(std::ifstream &in,
                            std::vector<std::vector<u8>> &values) {
  u64 size;
  readValue(in, size);
  values.resize(size);
  for (auto &value : values)
    readVector(in, value);
}

template <typename Imt>
void writeImt(std::ofstream &out, const Imt &imt) {
  writeValue(out, imt.mTableSize);
  writeValue(out, imt.mCmpsize);
  writeValue(out, imt.mDim);
  writeValue(out, imt.mDelta);
  writeValue(out, imt.mMetric);
  writeValue(out, imt.mCmp_len);
  writeValue(out, imt.mShareSize);
  writeValue(out, imt.mOutputSize);
  writeBitVector(out, imt.e);
  writeBitVector(out, imt.d);
  writeBitVector(out, imt.e_and);
  writeBitVector(out, imt.d_and);
  writeBitVector(out, imt.mImt_e_Share);
  writeBitVector(out, imt.mImt_d_Share);
}

template <typename Imt>
void readImt(std::ifstream &in, Imt &imt) {
  readValue(in, imt.mTableSize);
  readValue(in, imt.mCmpsize);
  readValue(in, imt.mDim);
  readValue(in, imt.mDelta);
  readValue(in, imt.mMetric);
  readValue(in, imt.mCmp_len);
  readValue(in, imt.mShareSize);
  readValue(in, imt.mOutputSize);
  readBitVector(in, imt.e);
  readBitVector(in, imt.d);
  readBitVector(in, imt.e_and);
  readBitVector(in, imt.d_and);
  readBitVector(in, imt.mImt_e_Share);
  readBitVector(in, imt.mImt_d_Share);
}

inline void writePeqt(std::ofstream &out, const PeqtSender &peqt) {
  writeValue(out, peqt.mDataSize);
  writeValue(out, peqt.mEqLength);
  writeBitVector(out, peqt.rShare);
  writeVector(out, peqt.tShare);
  writeNestedBytes(out, peqt.U);
  writeNestedBytes(out, peqt.V);
  writeVector(out, peqt.vose_val);
  writeBitVector(out, peqt.vose_table);
}

inline void readPeqt(std::ifstream &in, PeqtSender &peqt) {
  readValue(in, peqt.mDataSize);
  readValue(in, peqt.mEqLength);
  readBitVector(in, peqt.rShare);
  readVector(in, peqt.tShare);
  readNestedBytes(in, peqt.U);
  readNestedBytes(in, peqt.V);
  readVector(in, peqt.vose_val);
  readBitVector(in, peqt.vose_table);
}

inline void writePeqt(std::ofstream &out, const PeqtReceiver &peqt) {
  writeValue(out, peqt.mDataSize);
  writeValue(out, peqt.mEqLength);
  writeBitVector(out, peqt.rShare);
  writeVector(out, peqt.tShare);
  writeVector(out, peqt.eps);
  writeNestedBytes(out, peqt.W);
  writeBitVector(out, peqt.vose_table);
}

inline void readPeqt(std::ifstream &in, PeqtReceiver &peqt) {
  readValue(in, peqt.mDataSize);
  readValue(in, peqt.mEqLength);
  readBitVector(in, peqt.rShare);
  readVector(in, peqt.tShare);
  readVector(in, peqt.eps);
  readNestedBytes(in, peqt.W);
  readBitVector(in, peqt.vose_table);
}

inline std::filesystem::path freshPath(const std::string &root,
                                       const char *party) {
  return std::filesystem::path(root) / (std::string("l2_") + party + ".fresh");
}

inline void writeHeader(std::ofstream &out, const Key &key) {
  writeValue(out, kMagic);
  writeValue(out, kVersion);
  writeValue(out, key);
}

inline void readHeader(std::ifstream &in, const Key &expected) {
  u64 magic, version;
  Key actual{};
  readValue(in, magic);
  readValue(in, version);
  readValue(in, actual);
  if (magic != kMagic || version != kVersion ||
      std::memcmp(&actual, &expected, sizeof(Key)) != 0)
    throw std::runtime_error("L2 offline cache parameters do not match this run");
}

template <typename WriteFn>
void saveFresh(const std::filesystem::path &path, const Key &key,
               WriteFn &&writeFn) {
  std::filesystem::create_directories(path.parent_path());
  if (std::filesystem::exists(path))
    throw std::runtime_error("refusing to overwrite unconsumed L2 offline cache");
  const auto temporary = path.string() + ".tmp";
  std::ofstream out(temporary, std::ios::binary | std::ios::trunc);
  if (!out)
    throw std::runtime_error("cannot create L2 offline cache");
  writeHeader(out, key);
  writeFn(out);
  out.close();
  if (!out)
    throw std::runtime_error("failed to finalize L2 offline cache");
  std::filesystem::rename(temporary, path);
}

template <typename ReadFn>
void loadAndConsume(const std::filesystem::path &path, const Key &key,
                    ReadFn &&readFn) {
  if (!std::filesystem::exists(path))
    throw std::runtime_error("missing fresh L2 offline cache (it may already be consumed)");
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot open L2 offline cache");
  readHeader(in, key);
  readFn(in);
  in.close();
  if (!in)
    throw std::runtime_error("failed to finalize L2 offline cache read");
  std::filesystem::rename(path, path.string() + ".used");
}

} // namespace CmpFuzzyPSI::l2cache
