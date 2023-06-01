
require 'test/unit'
require 'fileutils'

TFILES = "tfiles/"
MAT1 = TFILES + 'file1.mata'
BINDIR = "../bin"

class MatConvertersTest < Test::Unit::TestCase
  def test_mata2mat
    tmpmata = TFILES + "trash.mata"
    tmpmat_out = TFILES + "trash.mat"
    FileUtils.cp MAT1, tmpmata
    pr = "mata2mat"
    system "#{BINDIR}/#{pr} #{tmpmata}"
    assert(File.exist?(tmpmat_out), "#{tmpmat_out} exists")
    arr = IO.read(tmpmat_out).unpack('iif*')
    exp = [4,3, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0]
    assert_equal(exp, arr)
    File.unlink(tmpmata)
    assert(!File.exist?(tmpmata), "#{tmpmata} does not exist")

    pr = "mat2mata"
    tmpmat = TFILES + 'trash.mat'
    tmpmat_out = TFILES + 'trash.mata'
    assert(File.exist?(tmpmat), "#{tmpmat} exists")
    system "#{BINDIR}/#{pr} #{tmpmat}"
    assert(File.exist?(tmpmat_out), "#{tmpmat_out} exists")
    arr = IO.readlines(tmpmat_out)
    assert_equal(%w(4 3), arr[0].split(" "))
    assert_equal(%w(1 2 3), arr[1].split(" "))
    assert_equal(%w(10 11 12), arr[4].split(" "))
      
    File.unlink(tmpmat)
    File.unlink(tmpmat_out)
  end

  
end



