
require 'test/unit'


$WIN32 = false

if ENV["OS"] =~ /Windows/
  $WIN32 = true
end

OBIWARP_PATH = "../bin/obiwarp" + ($WIN32 ? ".exe" : "")

TFILES = "tfiles/"
LMAT1 = TFILES + 'tmp1.lmat'
LMAT2 = TFILES + 'tmp1B.lmat'

class MyTests < Test::Unit::TestCase
  def diagnostics(reply)
    #puts reply
    hash = {}
    looking = false
    reply.split("\n").each do |line|
      if line =~ /\*{10,}/
        if looking
          looking = false
        else
          looking = true
        end
      elsif line =~ /(.*): (.*)/ && looking
        hash[$1.dup] = $2.dup
      end 
    end
    #p hash
    hash
  end
  
  def ob
    OBIWARP_PATH
  end

  def test_min_input
    #puts OBIWARP_PATH
    assert( File.exist?(ob), "obiwarp executable is in #{OBIWARP_PATH}")
    reply = `"#{ob}"`
    assert_match( /USAGE:/, reply, "no values passed in" )
    assert_match( /USAGE: #{File.basename(ob).gsub(/\.exe$/, '')}/, reply, "help progname matches executable")
    reply = `"#{ob}" only_1_file`
    assert_match( /USAGE:/, reply, "only one file passed in" )
  end

  def test_bad_files_input
    reply = `#{ob} badfile1 badfile2`
    assert_match(/Cannot open/, reply)
  end

  def test_opts
    t_opt(["--format"], "format", "mat");
    # What about no format give (should be same as LMAT1)
    t_opt_nil("format", "lmat")
    t_opt(%w(--outfile -o), "outfile", "myoutfilename");
    t_opt(%w(--images), "images");
    t_opt(%w(--timefile -t), "timefile", "mytimefile");

    expect = %w(cor cov prd euc)
    expect.each do |arg|
      t_opt(%w(--score -s), "score", arg)
    end
    t_opt(%w(--local -l), "local");
    t_opt(%w(--nostdnrm), "nostdnrm");
    t_opt_split(%w(--factor -f), "factor_diag", "factor_gap", "3.2,2.2")
    t_opt_split(%w(--gap -g), "gap_init", "gap_extend", "3.2,2.2")
    t_opt(%w(--init -i), "init_penalty", 2.1)
    t_opt(%w(--response -r), "response", 2.3);
  end

  ##########################################
  # HELPER FUNCS:
  ##########################################

  # For testing args like this: '23.3,5.2'
  def t_opt_split(opt_list, varname1, varname2, val)
    opt_list.each do |opt|
      reply =  `#{ob} #{opt} #{val} --diagnostics #{LMAT1} #{LMAT2}`
      hash = diagnostics(reply)
      val1, val2 = val.split(",")
      assert_equal("#{val1}", hash[varname1])
      assert_equal("#{val2}", hash[varname2])
    end
  end

  # for a variable we expect to see, even though no variables are passed in
  def t_opt_nil(varname, val)
    cmd = "#{ob} --diagnostics #{LMAT1} #{LMAT2}"
    reply = `#{cmd}`
    hash = diagnostics(reply)
    assert_equal("#{val}", hash[varname])
  end

  # for testing normal options
  # opt_list is a list of equivalent options
  # varname is the name of the diagnostic hash variable name
  # val is the value of the option passed in and the value expected out
  # if val == nil then the option is a flag and the output should be == 1
  def t_opt(opt_list, varname, val=nil)
    opt_list.each do |opt|
      cmd = "#{ob} #{opt} #{val} --diagnostics #{LMAT1} #{LMAT2}"
      #puts cmd
      reply = `#{cmd}`
      #puts "REPLY"
      #puts reply
      #puts "END PREPFYUDF"
      hash = diagnostics(reply)
      #p hash
      if val == nil
        assert_equal("1", hash[varname])
      else
        #puts "AHSHH" + "#{val}"
        #puts hash["score"]
        assert_equal("#{val}", hash[varname])
      end
    end
  end

end



