
require 'test/unit'
require 'ostruct'

$WIN32 = false; if ENV["OS"] == 'Windows_NT' then $WIN32 = true end

# updat_values copied from ruby facets:
# http://facets.rubyforge.org/doc/api/core/classes/Hash.html
class Hash
  def update_values
    each{ |k,v| store( k, yield(v) ) }
  end
end

OBIWARP_PATH = "../bin/obiwarp" + ($WIN32 ? ".exe" : '')
TFILES = "tfiles/"

hash = {
  :lmata1 => 'tmp1.lmata',
  :lmata1B => 'tmp1B.lmata',
  :lmata2 => 'tmp2.lmata',

  :lmat1 => 'tmp1.lmat',
  :lmat1B => 'tmp1B.lmat',
  :lmat2 => 'tmp2.lmat',
  :lmat_warped_default_G => 'tmp1B.lmat.warped_default',
  :lmat_warped_default => 'tmp1B.lmat.warped',

  :mat1 => 'tmp1.mat',
  :mat2 => 'file1.mat',
  :mat3 => 'file3.mat',
  :mat4 => 'file4.mat',
  :mata1 => 'tmp1.mata',
  :mata2 => 'file1.mata',
  :mata3 => 'file3.mata',
  :mata4 => 'file4.mata',
  :mat1_no_header => 'tmp1_no_header.mat',
  :mat2_no_header_messy => 'tmp1_no_header_messy.mat',
}.update_values {|v| TFILES + v }

F = OpenStruct.new(hash)

# Basic tests to ensure that files are being read, things are being warped
# when they are supposed to, etc.  Options are tested in cmdparser.
class ObiWarpTest < Test::Unit::TestCase

  @@lmat1_times = "1200.34 1212.34 1224.34 1236.34 1248.34 1260.34 1272.34 1284.34 1296.34 1308.34 1320.34 1332.34 1344.34 1356.34 1368.34 1380.34 1392.34 1404.34 1416.34 1428.34 1440.34 1452.34 1464.34 1476.34 1488.34 1500.34 1512.34 1524.34 1536.34 1548.34 1560.34 1572.34 1584.34 1596.34 1608.34 1620.34 1632.34 1644.34 1656.34 1668.34\n"
  @@mat1_times = "0 1 2 3 4 5\n"
  @@mat4_times = "0 1 2 3 4 5 6 7 8\n"
  def ob; OBIWARP_PATH end

  def test_min_input
    assert( File.exist?(ob), "obiwarp executable is in #{OBIWARP_PATH}")
    reply = `#{ob}`
    assert_match( /USAGE:/, reply, "no values passed in" )
    assert_match( /USAGE: #{File.basename(ob).gsub(/\.exe$/,'')}/, reply, "help progname matches executable")
    reply = `#{ob} only_1_file`
    assert_match( /USAGE:/, reply, "only one file passed in" )
  end

  # asserts that the file exists and is the same as "against" and deletes "file"

  def test_self_vs_self
    { F.mat1 => @@mat1_times,
      F.mata1 => @@mat1_times,
      F.lmat1 => @@lmat1_times,
      F.lmata1 => @@lmat1_times,
      F.mat4 => @@mat4_times,
    }.each do |k,v| 
      vs_self(k,v)
    end
  end

  def vs_self(file, expected)
    reply = `#{ob} #{file} #{file}`
    assert_equal(expected, reply)
  end

  def test_vs_other
    [
      [F.mat3, F.mat4, @@mat4_times],
      [F.mata3, F.mata4, @@mat4_times],
      [F.lmat2, F.lmat1, @@lmat1_times],
      [F.lmata2, F.lmata1, @@lmat1_times],
    ].each do |pair|
      assert_new_times(*pair)
    end
  end

  # returns true if the reply != the given times
  def assert_new_times(file1, file2, file2_times)
    reply = `#{ob} #{file1} #{file2}`
    #puts reply
    assert_equal(file2_times.chomp.split(" ").size, reply.chomp.split(" ").size, "same number of values")
    assert_not_equal(file2_times, reply, "times should not be the same after warping")
  end

  ##############################################
  ## HELPERS:
  ##############################################

end

